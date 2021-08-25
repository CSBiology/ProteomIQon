module Pipeline

open Expecto
open System
open System.IO

[<Tests>]
let pipelineTests =
    testList "Pipeline" [
        testCase "PeptideDB" <| fun _ ->
            let directory = Environment.CurrentDirectory
            printfn "%s" directory
            Expect.isTrue true "I compute, therefore I am."
    
        testCase "PeptideSpectrumMatching" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."

        testCase "PSMStatistics" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."

        testCase "PSMBasedQuantification" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."

        testCase "ProteinInference" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."
    ]
