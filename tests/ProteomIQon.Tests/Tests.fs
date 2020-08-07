namespace ProteomIQon.Tests

open Expecto
open FsCheck
open GeneratorsCode
open ProteomIQon
open ProteomIQon.TableSort
open ProteomIQon.Dto
open ProteomIQon.Domain
open System.IO

module Tests =
    let config10k = { FsCheckConfig.defaultConfig with maxTest = 10000}
    // bug somewhere:  registering arbitrary generators causes Expecto VS test adapter not to work
    //let config10k = { FsCheckConfig.defaultConfig with maxTest = 10000; arbitrary = [typeof<Generators>] }
    let configReplay = { FsCheckConfig.defaultConfig with maxTest = 10000 ; replay = Some <| (1940624926, 296296394) }

    [<Tests>]
    let testSimpleTests =

        testList "DomainTypes.Tag" [
            testCase "equality" <| fun () ->
                let result = 42
                Expect.isTrue (result = 42) "Expected True"

            testPropertyWithConfig config10k "whitespace" <|
                fun  () ->
                    Prop.forAll (Arb.fromGen <| whitespaceString())
                        (fun (x : string) -> 
                            x = x)
        ]

    [<Tests>]
    let tableSortTest =
        testList "DomainTypes.Tag" [
            testCase "equality" <| fun () ->
                let testfolder = (__SOURCE_DIRECTORY__ + @"/TableSortTest/")
                let quantFiles =
                    Directory.GetFiles(testfolder,("*.quant"))
                let protFiles =
                    Directory.GetFiles(testfolder,("*.prot"))
                let param = 
                    System.IO.File.ReadAllText(testfolder+"TableSortParams.json")
                    |> Json.deserialize<Dto.TableSortParams>
                    |> TableSortParams.toDomain
                sortTables quantFiles protFiles testfolder param
                let created = File.ReadAllLines (testfolder+"TableSort.tab")
                let reference = File.ReadAllLines (testfolder+"Reference.tab")
                let created_h = File.ReadAllLines (testfolder+"TableSort.tab")
                let reference_h = File.ReadAllLines (testfolder+"Reference.tab")
                Expect.isTrue (created = reference) "Expected True"
                Expect.isTrue (created_h = reference_h) "Expected True"
                File.Delete (testfolder+"TableSort.tab")
                File.Delete (testfolder+"TableSort_horizontal.tab")
                File.Delete (testfolder+"TableSort_log.txt")
                File.Delete (testfolder+"Test_log.txt")
        ]