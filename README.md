# Endonuclease-Digestion-Simulator-App
## Description

The **Endonuclease Digestion Simulator** is a user-friendly web application developed using Genie Builder and Julia. It is designed to help users simulate the digestion of DNA sequences with some of the most popular restriction enzymes. With this tool, users can input a DNA sequence and select up to three different restriction enzymes to visualize how the DNA is cut and the resulting fragments formed. This simulator is particularly valuable for molecular biologists, genetic engineers, researchers, and students interested in understanding and visualizing the effects of various restriction enzymes on DNA sequences.

### Key Features
- **Input DNA Sequences**: Easily enter your DNA sequence to be digested by the chosen enzymes.
- **Select Restriction Enzymes**: Pick from a variety of widely used restriction enzymes like EcoRI, HindIII, BamHI, AluI and more to model the cutting process.
- **Fragment Data Display**: Get detailed information on the resulting DNA fragments, including their sizes and the precise cut positions of each enzyme.
- **Simulated Gel Visualization**: The app generates a visual representation of the DNA fragments based on their sizes, as they would appear after gel electrophoresis.
- **Interactive Results**: View an interactive plot that shows the cut fragments and molecular weight markers, enabling a deeper exploration of the DNA digestion outcome.


https://github.com/user-attachments/assets/1e566aaf-3778-4821-b4ed-a095f7c6d477


### What are Restriction Enzymes?

Restriction enzymes, or **restriction endonucleases**, are specialized proteins that can perform cuts in DNA molecules at specific sites, known as recognition sites. These enzymes are naturally present in bacteria, serving as a defense mechanism against viral infections by interfering with viral replication chopping up that foreign DNA. In the field of molecular biology, restriction enzymes are essential for precision DNA manipulation, with key applications including:

- **Cloning**: Inserting a gene or a specific DNA segment into a vector for replication or expression in another organism.
- **DNA Mapping**: Studying the layout of restriction sites within a DNA sequence to gain insights into its structure or to compare it with other sequences.
- **Genetic Modification**: Making precise genetic changes in organisms by adding, removing, or replacing genes.

With the **Endonuclease Digestion Simulator**, users can explore the digestion process in an interactive manner, making it a highly effective tool for both education and research.

### What is Gel Electrophoresis?

**Gel electrophoresis** is a fundamental laboratory method used to separate nucleic acids or proteins by size and charge. For DNA analysis, this involves placing DNA fragments into a gel matrix and applying an electric field. Since DNA molecules are negatively charged, they migrate toward the positive electrode, with smaller fragments moving faster through the gel than larger ones. This produces a logarithmic pattern of bands that represents DNA fragments of varying sizes.

Gel electrophoresis is a essential technique in molecular biology for:

- **Analyzing DNA Fragments**: Visualizing and determining the sizes of DNA fragments produced from enzyme digestion.
- **Validating Cloning Procedures**: Checking whether a cloning attempt has successfully included the desired DNA fragment.
- **Forensic DNA Analysis**: Identifying individuals through their unique DNA profiles, as is done in forensic investigations.

The **Endonuclease Digestion Simulator** incorporates this concept by offering a simulated view of gel electrophoresis results, allowing users to observe how their digested DNA fragments would look on a gel, providing a comprehensive learning experience.

## Installation
Clone the repository and install the following dependencies:

+ GenieFramework
+ DataFrames
+ PlotlyBase
+ StipplePlotly

First cd into the project directory then run:
```
$> julia --project -e 'using Pkg; Pkg.instantiate()'

```
Finally, run the app

```
$> julia --project
```

```
julia> using GenieFramework
julia> Genie.loadapp() # load app
julia> up() # start server
```
## Usage

Open your browser and navigate to http://localhost:8000/.
