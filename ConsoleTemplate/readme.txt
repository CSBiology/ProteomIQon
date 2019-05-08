##install template
dotnet new --install C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\ConsoleTemplate\template

##verify installation --> look for pct 
dotnet new --list 

##create project from template in src folder, exchange "projectName" with whatever you want
C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src dotnet new pct -n "projectName" --force

## add project to solution
C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon>dotnet sln ProteomIQon.sln add "C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src\test3\test3.fsproj"