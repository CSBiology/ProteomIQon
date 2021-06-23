(**
// can't yet format YamlFrontmatter (["title: How to document your work"; "category: Developer Notes"; "categoryindex: 2"; "index: 2"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 8 }) to pynb markdown

# How to document your work
**Disclaimer** This page gives only an overview of the FSharp.Formatting functionality needed to contribute to the ProteomIQon, for details please visit 
the documentation of [FSharp.Formatting](https://fsprojects.github.io/FSharp.Formatting/).

The web page that you are just reading was created using the [fsdocs](https://fsprojects.github.io/FSharp.Formatting/commandline.html) dotnet tool!
Thanks to the effort of the FSharp.Formatting team documenting the ProteomIQon is a straightforward process.

## Building the docs:
Assuming that you have cloned the ProteomIQon repository you can build the docs by navigating to the project root and calling: 

	dotnet fsdocs watch --eval --noapidocs

Executing this command will lead to parsing of all .fsx files in the docs folder, this will finally lead to formatting of every .fsx as a html documents and notebooks.
This is all done behind the scenes orchestrated by the fsdocs cli tool. Once this process is finished your default browser should start and automatically navigate to 
the adress of the index.html hosted on a local webserver.

## How to change exisiting content:
The fsdocs tool enables you to manipulate the docs and immediately observe the effects of your actions. To try out the hot reload functionality simply change the content of a 
fsx file placed in the ./docs folder of the proteomiqon repository and save your changes. After a short while, this should trigger a reload in you webbrowser and you should be able to 
see an updated .html file with your changes incorporated.

## How to add a tool documentation:
If you are working on your own tool, chances are high that you want to add a new .html to the docs. A suggested workflow to add a docs page for your tool could look like this:

1. Navigate to ./docs/tools
2. Copy one of the existing .fsx files and rename it
3. Replace the content with the content fitting to your tool.

## Update the online documentation:
Once you are happy with your contribution simply open a pull request. Upon acceptance of your commit, a github action is triggered and the new version of the docs are automatically released.

*)