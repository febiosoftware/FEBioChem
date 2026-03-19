## From FEBio Studio
You can download and install the FEBioChem plugin in FEBio Studio using the Plugin Repo. 

## From the command line
Like any other FEBio plugin, the plugin must be placed in a folder and the path to the plugin must be defined in the FEBio configuration file. This file is usually called `febio.xml` and can be found in the same location as the FEBio executable. In this file, add the following line:

```xml
<import>C:\path\to\febio\plugin\FEBioChem.dll</import>
```

Make sure to include the full path name and file name of the plugin. When FEBio starts it will read the configuration file and load all the plugins defined therein. A message will be shown to the screen to inform the user whether the plugin was loaded successfully or not.