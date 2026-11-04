namespace ProteomIQon

/// Processing of the raw XIC of a precursor.
module Xic =

    /// A zero point that follows a point more than 100 times higher is a drop-out; the mask is false there.
    let keepMask (y: float []) =
        Array.init y.Length (fun i ->
            if i = 0 || i = y.Length - 1 || y.[i] > 0. then
                true
            else
                let previous = y.[i - 1]
                previous = 0. || previous <= 100. * (y.[i] + 1.)
            )

    /// ALS baseline subtraction for XICs of 3 to 500 points; negative values are set to 0.
    let subtractBaseline (baseline: Domain.BaseLineCorrection option) (y: float []) =
        match baseline with
        | Some p when y.Length >= 3 && y.Length <= 500 ->
            let b = BaselineAls.baselineAls p.MaxIterations p.Lambda p.P y
            Array.map2 (fun v bv -> max 0. (v - bv)) y b
        | _ -> Array.copy y
