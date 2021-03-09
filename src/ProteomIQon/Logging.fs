namespace ProteomIQon

open NLog
open NLog.Targets
open NLog.Config
open System.IO

module Logging =

    let generateConfig (folderPath: string)= 
        //initializes base configuration class, can be modified
        let config = new LoggingConfiguration()

        //initialises base console target, can be modified
        let consoleTarget = new ColoredConsoleTarget("console")
        //new parameters for console target
        let layoutConsole = new Layouts.SimpleLayout (@"${date:format=HH\:mm\:ss} ${level} ${message} ${exception}")
        consoleTarget.Layout <- layoutConsole

        //initialises base file target, can be modified
        let fileTarget = new FileTarget("file")
        //new parameters for file target
        let fileName = new Layouts.SimpleLayout (Path.Combine (folderPath, @"${logger}_log.txt"))
        let layoutFile = new Layouts.SimpleLayout ("${longdate} ${logger} ${level:uppercase=true} ${message} ${exception}")
        fileTarget.FileName <- fileName
        fileTarget.Layout <- layoutFile

        config.AddTarget(consoleTarget)
        config.AddTarget(fileTarget)

        //declares which error results in a log in which target
        config.AddRuleForOneLevel(LogLevel.Trace, fileTarget)
        config.AddRuleForOneLevel(LogLevel.Info, consoleTarget)

   
        //activates config for logger
        LogManager.Configuration <- config

    let createLogger (loggerName: string) = 


        //new instance of "Logger" with activated config
        let logger = LogManager.GetLogger(loggerName)

        logger