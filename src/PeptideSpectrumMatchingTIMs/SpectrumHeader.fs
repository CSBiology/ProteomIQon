namespace ProteomIQon

// A streaming parser of the spectrum description that reads only the fields the search needs.
module TimSpectrumHeader =

    open Newtonsoft.Json

    type SpectrumHeader =
        {
            ID          : string
            /// MS:1000511 of the spectrum, None when absent.
            MsLevel     : int option
            /// MS:1000016 of the first scan that carries it.
            ScanTime    : float option
            /// Of the first selected ion that carries one of the two.
            PrecursorMz : float option
            /// MS:1000041 of the first selected ion that carries it.
            ChargeState : int option
            /// MS:1002815 of the first scan that carries it.
            IonMobility : float option
            /// The message of the failure when the description could not be parsed, else null.
            /// Such a header carries no other value.
            Error       : string
        }

    let private parseFloat (text: string) =
        match System.Double.TryParse(text, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture) with
        | true, v -> Some v
        | _ -> None

    /// Reads the first entry of the Values array of the CV parameter object the reader stands in,
    /// and leaves the reader at the end of that object.
    let private firstValue (reader: JsonTextReader) =
        let mutable value : string = null
        let mutable depth = 0
        let mutable started = false
        while (not started || depth > 0) && reader.Read() do
            if reader.TokenType = JsonToken.StartObject then started <- true
            match reader.TokenType with
            | JsonToken.StartObject -> depth <- depth + 1
            | JsonToken.EndObject -> depth <- depth - 1
            | JsonToken.PropertyName when depth = 1 && (reader.Value :?> string) = "Values" ->
                if reader.Read() && reader.TokenType = JsonToken.StartArray then
                    if reader.Read() && reader.TokenType <> JsonToken.EndArray then
                        value <- (match reader.Value with null -> null | v -> System.Convert.ToString(v, System.Globalization.CultureInfo.InvariantCulture))
                        // skip the rest of the array
                        let mutable arrayDepth = 1
                        while arrayDepth > 0 && reader.Read() do
                            match reader.TokenType with
                            | JsonToken.StartArray -> arrayDepth <- arrayDepth + 1
                            | JsonToken.EndArray -> arrayDepth <- arrayDepth - 1
                            | _ -> ()
            | _ -> ()
        value

    /// Parses one description. The path of property names from the root tells where a CV
    /// parameter sits. A parameter of the spectrum itself has the shortest path. A scan and a
    /// selected ion each add their own levels to it.
    let parse (json: string) =
        use reader = new JsonTextReader(new System.IO.StringReader(json))
        reader.DateParseHandling <- DateParseHandling.None
        reader.FloatParseHandling <- FloatParseHandling.Double
        let path = System.Collections.Generic.List<string>()
        let mutable rootStarted = false
        let mutable rootClosed = false
        let mutable pending : string = null
        let mutable id : string = null
        let mutable msLevel = None
        let mutable scanTime = None
        let mutable mobility = None
        let mutable precursorMz = None
        let mutable charge = None
        // per selected ion
        let mutable ionPrecursorMz = None
        let mutable ionSelectedMz = None
        // path.[0] is the root object
        let inSelectedIon () =
            path.Count = 8 && path.[1] = "Precursors" && path.[2] = "properties" && path.[4] = "SelectedIons" && path.[5] = "properties" && path.[7] = "properties"
        let inScan () =
            path.Count = 5 && path.[1] = "Scans" && path.[2] = "properties" && path.[4] = "properties"
        let inSpectrum () = path.Count = 2 && path.[1] = "properties"
        let closeIon () =
            if precursorMz.IsNone then
                precursorMz <- (match ionPrecursorMz with Some v -> Some v | None -> ionSelectedMz)
            ionPrecursorMz <- None
            ionSelectedMz <- None
        while reader.Read() do
            match reader.TokenType with
            | JsonToken.PropertyName ->
                let name = reader.Value :?> string
                pending <- name
                if path.Count = 1 && name = "ID" then
                    if reader.Read() then id <- (reader.Value :?> string)
                elif inSpectrum () && name = "MS:1000511" then
                    let v = firstValue reader
                    if msLevel.IsNone && not (isNull v) then msLevel <- parseFloat v |> Option.map int
                elif inScan () && name = "MS:1000016" then
                    let v = firstValue reader
                    if scanTime.IsNone && not (isNull v) then scanTime <- parseFloat v
                elif inScan () && name = "MS:1002815" then
                    let v = firstValue reader
                    if mobility.IsNone && not (isNull v) then mobility <- parseFloat v
                elif inSelectedIon () && name = "MS:1002234" then
                    let v = firstValue reader
                    if not (isNull v) then ionPrecursorMz <- parseFloat v
                elif inSelectedIon () && name = "MS:1000744" then
                    let v = firstValue reader
                    if not (isNull v) then ionSelectedMz <- parseFloat v
                elif inSelectedIon () && name = "MS:1000041" then
                    let v = firstValue reader
                    if charge.IsNone && not (isNull v) then charge <- parseFloat v |> Option.map int
            | JsonToken.StartObject ->
                if path.Count = 0 then
                    if rootStarted then failwith "More than one JSON root object."
                    rootStarted <- true
                path.Add(if isNull pending then "" else pending)
                pending <- null
            | JsonToken.EndObject ->
                if inSelectedIon () then closeIon ()
                if path.Count > 0 then path.RemoveAt(path.Count - 1)
                if path.Count = 0 then rootClosed <- true
            | JsonToken.StartArray ->
                // arrays outside Values are not expected; skip them whole
                let mutable depth = 1
                while depth > 0 && reader.Read() do
                    match reader.TokenType with
                    | JsonToken.StartArray -> depth <- depth + 1
                    | JsonToken.EndArray -> depth <- depth - 1
                    | _ -> ()
                pending <- null
            | _ -> pending <- null
        if not rootStarted || not rootClosed || path.Count <> 0 then
            failwith "Incomplete spectrum description JSON."
        {
            ID = id
            MsLevel = msLevel
            ScanTime = scanTime
            PrecursorMz = precursorMz
            ChargeState = charge
            IonMobility = mobility
            Error = null
        }

    /// parse, with a failure turned into a header whose Error names it.
    let tryParse (json: string) =
        try parse json
        with ex -> { ID = null; MsLevel = None; ScanTime = None; PrecursorMz = None; ChargeState = None; IonMobility = None; Error = ex.Message }

