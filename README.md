# FBXToBoneMF

Converts FBX Files to the Bonetown Model Format for use in Minecraft

## Usage:

./FBXToBoneMF.exe [fbx file] [name of output file] [json | cbor]

There are 2 options for export, cbor is the one used by the Bonetown mod itself. Json is a human readable version of the cbor 
so you can inspect the output.

Example usage:
`
./FBXToBoneMF.exe biped.fbx biped cbor
`

## Current Runs On:

Windows

## Dependencies:

* [FBX SDK](https://www.autodesk.com/developer-network/platform-technologies/fbx-sdk-2020-0)
