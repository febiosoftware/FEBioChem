## Using the FEBioChem plugin
Like any other FEBio plugin, the plugin must be placed in a folder and the path to the plugin must be defined in the FEBio configuration file. This file is usually called `febio.xml` and can be found in the same location as the FEBio executable. In this file, add the following line:

```xml
<import>C:\path\to\febio\plugin\FEBioChem.dll</import>
```

Make sure to include the full path name and file name of the plugin. When FEBio starts it will read the configuration file and load all the plugins defined therein. A message will be shown to the screen to inform the user whether the plugin was loaded successfully or not.

## Defining an FEBioChem model
For the most part an FEBioChem model is a standard FEBio model with a few modifications. These are discussed in the following sections. Note that this user’s manual assumes version 4.0 of the FEBio file format.

### Module
The `type` attribute of `Module` tag needs the value `reaction-diffusion`.

```xml
<Module type="reaction-diffusion"/>
```

This will inform FEBio to load the reaction-diffusion solver that is part of the FEBioChem plugin.

This module does not include convection. If you want to include the convection term, you need to use the `reaction-diffusion-convection` module.

```xml
<Module type="reaction-diffusion-convection"/>
```

Now, a velocity field variable is added to the model, which can be defined in the `Initial` section of the input file. See section [initial velocity](#initial-velocity) for more information on how to initialize the velocity field.

### Globals
The chemical species that are used in the model are defined in the `Globals` section. The species are defined identically to the way they are defined for a multiphasic problem and we refer to the FEBio user’s manual. The following caveat applies:

* A **solute** defines a chemical species that is not bound to the solid phase of a mixture. If a model doesn’t represent a mixture, each solute defines a chemical species. No child properties need to be defined for a solute.
* A **solid_bound** defines a chemical species that is bound to the solid-phase of a mixture. Again, no child properties need to be defined for a solid bound species.

## Solver
Make sure to set the solver type to the same name as the module. For instance, for a reaction-diffusion model, set

```xml
<solver type="reaction-diffusion">
...
</solver>
```

### Material
The material type needed by this plugin is `reaction-diffusion`. It needs the following material properties. Note that not all are required.

|Property|Description|
|--------|-----------|
|`solid_volume_fraction` | The volume fraction of the solid phase|
|`species` | Specifies the properties of a chemical species|
|`solid_bound_species` | Specifies a species that remains bound to the solid phase of the mixture.
|`reaction`| Specifies the chemical reaction.|

The `solid_volume_fraction` is only required when modeling a mixture. It has an initial value of zero.

#### species property
Each chemical species that is used by the material must be defined in the material’s definition. (Note that this can be less than the total number of species in the model, as defined by the `Globals` section.). A particular species is referenced by name. For each species, the `diffusivity` parameter must be defined.

```xml
<species name="A">
    <diffusivity>0.05</diffusivity>
</species>
```

#### solid_bound_species property
Each solid-bound chemical species that is used by the material must be included in the material’s definition. (Note that this can be less than the total number of species in the model, as defined by the `Globals` section). For each solid-bound species, define the following material parameters.

|Parameter | Description | Default|
|----------|-------------|--------|
|`rho0` |	Initial apparent density (w.r.t. the mixture volume) |	0|
|`density`| The true density of the solid-bound species | 0 |
|`rhomin` | The minimum value of the apparent density allowed	| 0 |
|`rhomax` | The maximum value of the apparent density allowed |	(ignored) |
|`molar_mass`| The molar mass of the solid-bound species. |0 |

#### reaction property

The `reaction` property allows you to define a chemical reaction. Note that currently FEBioChem only models forward reactions. To model a reversible reaction you must specify two reactions, where the first reaction defines the forward reaction and the second defines the reverse reaction. Defining a reaction in the FEBioChem plugin is very easy since only two parameters need to be specified.

The reaction equation is defined by writing the reaction formula as follows.

* Use the name of the species as defined in the Globals section to reference it.
* A stoichiometric coefficient of 1 does not need to be defined, but otherwise it must an integer number that appears before the species’ name and separated from it via a multiplication symbol `*`. 
* A plus symbol `+` is used to enumerate all reactants and products.
* The symbols `-->` or `->` can be used to separate the reactants from the products. 
* The reactant or product lists can be empty. 

Here are some examples. Assume that all the chemical species are properly defined in the `Globals` section. 

```
A-->B
A+B-->C
2*A+2*B->C+B
A-->
-->B
```

Again, it is emphasized that all species involved in a reaction must be defined inside the material’s definition.

### Initial values
For each species (but not solid-bound species), an initial value can be set in the `Initial` section of the FEBio input file via the `init` tag.

```xml
<ic type="initial concentration" node_set="set1">
	<dof>c1</dof>
	<value>1.0</value>
</init>
```

The `dof` parameter references the particular species (or solute). In this case, `c1` references the solute with ID equal to 1. The `node_set` references a node set that needs to be defined in the `Mesh` section of the input file.

### Boundary conditions
The syntax for defining fixed or prescribed is identical for any FEBio model. Use `c[n]` to reference a particular species (or solute, again no solid-bound species). Thus, for example, `c1` references solute with ID one, `c2` is the solute with ID 2, and so forth. 
An example of a fixed boundary condition looks like this.

```xml
<bc type="zero concentration" node_set="FixedSet1">
	<dof>c1</dof>
</bc>
```
An example of a prescribed boundary condition looks like this.

```xml
<bc type="prescribed concentration" node_set="PrescribedSet1">
	<dof>c1</dof>
	<value lc="1">1.0</value>
</prescribe>
```

The `node_set` attribute references nodesets that need to be defined in the `Mesh` section of the input file.

### Prescribing a concentration flux
A concentration flux can be prescribed via a surface load of type concentration flux and is defined as follows.

```xml
<surface_load type="concentration flux" surface="SurfaceLoad1">
	<solute_id>1</solute_id>
	<flux lc="1">-1.0</flux>
</surface_load>
```

The `surface` attribute references a surface that is defined in the `Mesh` section. The `solute_id` references the ID of the corresponding solute as defined in the `Globals` section of the input file.

### Plot Variables
The FEBioChem plugin defines several new plot variables. They are list in the table below.

|Plot variable |Description|
|--------------|-----------|
|`concentration` / `actual concentration` |	Get the true concentration of a species |
|`effective concentration` | Get the effective concentration of a species |
|`sbs concentration` | The true concentration of a solid-bound species|
|`sbs apparent density` | The apparent density of a solid-bound species|
|`solid volume fraction`| The solid volume fraction of the mixture |

Note that the variables that references a species or a solid-bound species need to specify the name of the corresponding species as part of the plot variable definition. For instance,

```xml
<var type="concentration['c1']"/>
```
As with any other plot variable definition in FEBio, you can define an alias.

```xml
<var type="concentration['c1']=c1/>
```

The alias is the name of the data field that will be stored in the plot file, and thus will be show in FEBioStudio.

### Initial velocity
For reaction-diffusion-convection problems, the fluid velocity needs to be initialized. This is done in the `Initial` section of the FEBio input file.

The following example sets the x-component of the velocity vector to 1.0.

```xml
<ic type="initial velocity" node_set="SomeNodeSet">
	<dof>vx</dof>
	<value>1.0</value>
</init>
```

You can also define a heterogeneous field. For instance,

```xml
<ic type="initial velocity" node_set="SomeNodeSet">
	<dof>vx</dof>
	<value type="map">init_vel<value/>
</init>
```
Here, `init_vel` is a node data field that is defined in the `MeshData` section.

You can also use a mathematical expression to set the initial velocity field. 

```xml
<ic type="initial velocity" node_set="SomeNodeSet">
	<dof>vx</dof>
	<value type="math">4*(0.5 - Y)*(Y + 0.5)<value/>
</init>
```
