namespace ProteomIQon

open System
open MzIO
open FSharpAux
open ProteomIQon.Domain
open MzTABSections
open MzTABFormater

module MzTAB =

    let createMzTab path tab prot quant qpsm (param: Domain.MzTABParams) =
        let protFiles =
            IO.Directory.GetFiles(prot,("*.prot"))
        let quantFiles =
            IO.Directory.GetFiles(quant,("*.quant"))
        let qpsmFiles =
            IO.Directory.GetFiles(qpsm,("*.qpsm"))
        let allAligned =
            allignAllFiles tab protFiles quantFiles qpsmFiles
        let protSection =
            proteinSection allAligned param
        let pepSection =
            peptideSection allAligned param
        let psmSection' =
            psmSection allAligned param
        metaDataSection path param.MetaData
        protHeader path param
        protBody path protSection
        pepHeader path param
        pepBody path pepSection
        psmHeader path param
        psmBody path psmSection'