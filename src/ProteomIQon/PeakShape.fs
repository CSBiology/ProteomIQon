namespace ProteomIQon

open System
open FSharp.Stats

/// Position and width of a fitted chromatographic peak. A peak is a Gaussian (parameters amplitude, mean, sigma) or an
/// exponentially modified Gaussian (amplitude, mean, sigma, tau); for the latter mean and sigma belong to the Gaussian before the
/// exponential tail, so the apex lies right of mean and the peak is wider than sigma.
module PeakShape =

    let private isUsable (v: float) = not (Double.IsNaN v || Double.IsInfinity v)

    /// The scaled complementary error function exp(t^2) * erfc(t). The error function of FSharp.Stats has an absolute error of about
    /// 1e-7, too much for the tail, so it is computed here: for t < 0 by erfcx(t) = 2 exp(t^2) - erfcx(-t), for 0 <= t < 2 as
    /// exp(t^2) - 2 / sqrt pi * (t + 2 t^3 / 3 + 4 t^5 / 15 + ...), from there on by the continued fraction
    /// 1 / sqrt pi / (t + (1/2) / (t + (2/2) / (t + (3/2) / (t + ...)))), which does not overflow for large t.
    let rec erfcx (t: float) =
        if Double.IsNaN t then
            nan
        elif t < 0. then
            2. * exp (t * t) - erfcx (-t)
        elif t < 2. then
            let t2 = t * t
            let mutable term = t
            let mutable sum = t
            let mutable n = 0
            while term > 1e-18 * sum && n < 200 do
                term <- term * 2. * t2 / float (2 * n + 3)
                sum <- sum + term
                n <- n + 1
            exp t2 - 2. / sqrt Math.PI * sum
        else
            let mutable f = t
            for k in 80 .. -1 .. 1 do
                f <- t + 0.5 * float k / f
            1. / (sqrt Math.PI * f)

    /// (mode - mean) / sigma of an exponentially modified Gaussian with r = tau / sigma. The mode solves
    /// erfcx(t) = r * sqrt(2 / pi) with t = (1 / r - u) / sqrt 2 and u = (mode - mean) / sigma. erfcx decreases in t, so the root is
    /// bracketed and found by Newton steps (erfcx'(t) = 2 t erfcx(t) - 2 / sqrt pi); a step that is not finite or leaves the bracket
    /// is replaced by a bisection step. 0 for a non-finite or non-positive r.
    let emgModeOffset (r: float) =
        if not (isUsable r && r > 0.) then
            0.
        else
            let target = r * sqrt (2. / Math.PI)
            // erfcx(t) < 1 / (t sqrt pi) for t > 0, so erfcx(hi) < target; erfcx(-t) >= exp(t^2) for t > 0 and erfcx(-1) > 1, so erfcx(lo) > target
            let mutable hi = 1. / (target * sqrt Math.PI) + 1.
            let mutable lo = -(sqrt (log (max target 1.)) + 1.)
            // for a small target erfcx(t) is close to 1 / (t sqrt pi)
            let mutable t = if target < 0.2 then 1. / (target * sqrt Math.PI) else 0.5 * (lo + hi)
            let mutable running = true
            let mutable steps = 0
            while running && steps < 200 do
                let value = erfcx t
                if value > target then lo <- t else hi <- t
                let newton = t - (value - target) / (2. * t * value - 2. / sqrt Math.PI)
                let next = if isUsable newton && newton > lo && newton < hi then newton else 0.5 * (lo + hi)
                running <- abs (next - t) >= 1e-15 * max 1. (abs t) && next <> lo && next <> hi
                t <- next
                steps <- steps + 1
            1. / r - sqrt 2. * t

    /// The model has an exponential tail if it has at least 4 parameters with a finite sigma > 0 and tau > 0.
    let private hasTail (peakParams: float []) =
        peakParams.Length >= 4 && isUsable peakParams.[2] && peakParams.[2] > 0. && isUsable peakParams.[3] && peakParams.[3] > 0.

    /// Scan time of the highest point of the fitted model: for an EMG mean + sigma * emgModeOffset (tau / sigma), for a Gaussian
    /// (3 parameters, or 4 with unusable tau) its mean, nan with fewer parameters.
    let apex (peakParams: float []) =
        if peakParams.Length < 3 then nan
        elif hasTail peakParams then peakParams.[1] + peakParams.[2] * emgModeOffset (peakParams.[3] / peakParams.[2])
        else peakParams.[1]

    /// Standard deviation of the fitted model: for an EMG sqrt (sigma^2 + tau^2), for a Gaussian (3 parameters, or 4 with unusable tau)
    /// sigma, nan with fewer parameters.
    let width (peakParams: float []) =
        if peakParams.Length < 3 then nan
        elif hasTail peakParams then sqrt (peakParams.[2] * peakParams.[2] + peakParams.[3] * peakParams.[3])
        else peakParams.[2]
