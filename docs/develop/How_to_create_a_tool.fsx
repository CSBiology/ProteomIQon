(**
---
title: How to your own tool
category: Developer Notes
categoryindex: 2
index: 1
---
*)
(**
# How to create a tool

If you want to extend the functionality of the ProteomIQon you can either submit a pull request to update
an existing tool or - if the functionality is better encapsulated as a stand alone application - create a new tool.
To make this as as easy as possible we created a tool template which sets up a basic scaffold for your new tool!

## Installing the template:
Assuming that you have cloned the ProteomIQon repository you can install the tool by navigating to the project root and calling: 
*)

(**
	dotnet new --install ./ConsoleTemplate/template
*)

(**
The installation can then be verified by executing the following snippet:
*)

(**
	dotnet new --list 
*)

(**
This should print a list of installed templates. If everything was successful you will find the line among the installed tools: 
*)

(**
	Templates                           Short Name               Language          Tags
	------------------------------      -------------------      ------------      ----------------------------------------
	...
	prototypical proteomiqon co...      pct                      F#                proteomiqon console/proteomiqon/template
	...
*)

(**
## Adding a new tool to the project:
To initialize a new tool using the template, navigate into the 'src' folder, replace 'projectName' with your choice for a tool name and call:
*)

(**
	dotnet new pct -n "projectName" --force
*)

(**
Afterwards navigate back to the project root and execute the following line to add your tool to the ProteomIQon solution:
*)

(**
	dotnet sln ProteomIQon.sln add "./src/projectName/projectName.fsproj"
*)

(**
Afterwards you should be good to go! Have fun extending the ProteomIQon, we look forward to your contribution!
*)
