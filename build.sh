#!/usr/bin/env bash
clear

dotnet run --project ./build/build.fsproj "$@"
