# OrbAlign
**Orbital Alignment**

A python-based quantum chemistry bonding analysis tool to align the occupied orbitals of a molecule towards its fragments so that trivial orbitals (that originates from isolated fragments, e.g. core electrons and inert bonds) will be labelled as "maximially overlapped orbitals", while non-trivial orbitals (that originates from fragment interactions or involves additional electrons) will be labelled as "minimally overlapped orbitals".

## Requirement
- Python 3
- Numpy 1.18.5
- Gaussian 09 or 16

This combination has been well tested. Other versions are not guaranteed to work but are welcomed to test.

## Tutorial
The overall workflow of OrbAlign basically involves the following steps:
- Perform quantum chemistry calculation on the molecule of interest
- Divide the molecule into fragments
- Perform quantum chemistry calculation on each fragment
- Align

In the following tutorial, I will use the cluster compound [Au4(PH3)4]2+ as an example. All mentioned files can be found in the `test` directory under the same root of this README file.

1. Perform Gaussian calculation on the molecule of interest

    First, optimize the molecule and prepare a Gaussian input file for single-point calculation with "nosymm" keyword so that the orientation of the molecule (and fragments in later steps) will be frozen. The keyword "pop=nboread", along with a seperated line `$nbo archive filename=Au4L4-nosymm $end` in the end of file, is also necessary for Gaussian to prepare a .47 file that contains relevant information. Run Gaussian and prepare the corresponding .fchk file.

    **Input**
    - *Au4L4-nosymm.gjf*
    ```
    %chk=Au4L4-nosymm.chk
    #p pbe1pbe/def2TZVP pop=nboread nosymm

    Au4L4-nosymm

    2 1
     Au                 0.98780500    0.98780500    0.98780500
     Au                -0.98780500   -0.98780500    0.98780500
     Au                 0.98780500   -0.98780500   -0.98780500
     Au                -0.98780500    0.98780500   -0.98780500
     P                  2.33876700    2.33876700    2.33876700
     H                  3.21542500    1.65931900    3.21542500
     H                  3.21542500    3.21542500    1.65931900
     H                  1.65931900    3.21542500    3.21542500
     P                 -2.33876700   -2.33876700    2.33876700
     H                 -3.21542500   -3.21542500    1.65931900
     H                 -1.65931900   -3.21542500    3.21542500
     H                 -3.21542500   -1.65931900    3.21542500
     P                  2.33876700   -2.33876700   -2.33876700
     H                  3.21542500   -3.21542500   -1.65931900
     H                  1.65931900   -3.21542500   -3.21542500
     H                  3.21542500   -1.65931900   -3.21542500
     P                 -2.33876700    2.33876700   -2.33876700
     H                 -1.65931900    3.21542500   -3.21542500
     H                 -3.21542500    1.65931900   -3.21542500
     H                 -3.21542500    3.21542500   -1.65931900

     $nbo archive filename=Au4L4-nosymm $end


    ```

    **Command**
    ```
    g16 Au4L4-nosymm.gjf Au4L4-nosymm.log
    formchk Au4L4-nosymm.chk Au4L4-nosymm.fchk
    ```

    **Output**
    - *Au4L4-nosymm.log*
    - *Au4L4-nosymm.fchk*
    - *Au4L4-nosymm.47*

    Depending on the NBO version you use in Gaussian, the output filename might be in upper case, i.e. `AU4L4-NOSYMM.47`, in which case you should manually rename it into lower case. Same for all Gaussian calculations below.

2. Divide the molecule into fragments

    Then, divide the molecule into several chemical fragments. Each fragment should be chemically meaningful.
    For example, in the [Au4(PH3)3]2+ cluster, we would like to know how the four Au atoms bond to each other. Apparently, their bonding is not a typical 2-center-2-electron bond. Hence we take each [Au]+ as a fragment; each [PH3] ligand is also a chemically meaningful entity and we expect they retain their electronic structures when forming the cluster (because they are "ligands": we expect each ligand donate one pair of electrons to an Au ion, but this electron pair should still be associated with the ligand, with only minimal deformation during the interaction with Au).
    Based on the fragmentation proposed above, we prepare a file recording the atom ids of each fragment, as well as its charge, multipliticy and fragment labelling.

    **Input**
    - *Au4L4.id*
    ```
    # Atom Ids  # Charge    # Multiplicity  # Fragment label    # Fragment .47 filename
    1           1           1               Au                  AU.47
    2           1           1               Au                  AU.47
    3           1           1               Au                  AU.47
    4           1           1               Au                  AU.47
    5-8         0           1               L
    9-12        0           1               L
    13-16       0           1               L
    17-20       0           1               L
    ```
    In this .id file, the first column records the atom ids of each fragment according to the atom numbering of the Gaussian input file `Au4L4-nosymm.gjf`.

    The second and third column record the charge and multiplicity of each fragment. In this case, Au atoms have +1 charge, PH3 ligands are neutral, both adopting singlet states.

    The 4th column is a label for each fragment. Identical labels (e.g. "L") are allowed and will be converted to numbered labels (e.g. "L-1", "L-2", "L-3" and "L-4") automatically if the fragment filenames are not explicitly provided.

    The 5th column denotes the .47 filename of each fragment. This column can always be left blanck safely. It's only useful when a fragment contains only a single atom, in which case one can perform a single run for this type of atom and replicates the result for all identical atoms.

    Each row represents a fragment unless it startswith a "#" in which case this row will be treated as omment line and be ignored.

    **Command**
    ```
    python autofrag.py Au4L4-nosymm.gjf Au4L4.id
    ```
    This script will automatically generate the Gaussian input files for each fragment except those with pre-specified .47 filenames (all Au fragments in our example). As mentioned above, you can always ignore this usage and leave blanck the last column of Au4L4.id, in which case there will be four more .gjf files generated in this step.

    A `Au4L4.id.rev` file will also be generated as an updated .id file, which have all .47 filenames explicitly specified now. By default, the filenames are given according to the fragment labels and charges. You don't have to understand or modify the content of this file.

    **Output**
    - *Au4L4.id.rev*
    - *Au4L4-nosymm-L-1-p0.gjf*
    - *Au4L4-nosymm-L-2-p0.gjf*
    - *Au4L4-nosymm-L-3-p0.gjf*
    - *Au4L4-nosymm-L-4-p0.gjf*

3. Perform quantum chemistry calculation on each fragment

    Now, run Gaussian calulation on each fragment. The Gaussian input files automatically prepared in the last step should just work, but it's still recommended to take a look at their contents to make sure they are legal Gaussian input files.

    **Input**
    - *Au-d10.gjf*
    - *Au4L4-nosymm-L-1-p0.gjf*
    - *Au4L4-nosymm-L-2-p0.gjf*
    - *Au4L4-nosymm-L-3-p0.gjf*
    - *Au4L4-nosymm-L-4-p0.gjf*

    **Command**
    ```
    for f in Au4L4-nosymm-*gjf;  do g16 $f; done
    ```

    **Output**
    - *Au-d10.47*
    - *Au4L4-nosymm-L-1-p0.47*
    - *Au4L4-nosymm-L-2-p0.47*
    - *Au4L4-nosymm-L-3-p0.47*
    - *Au4L4-nosymm-L-4-p0.47*

    Among all output files from Gaussian, only the .47 files are useful in the following step. If you prefer your own .47 filenames when performing Gaussian calculations, change them back to these default names or make sure they are consistent with the filenames in the `Au4L4.id.rev` file so that the program can still properly recognize them.
    It should be noted that the `Au-d10.gjf` job contains extra keywords to ensure the obtained wavefunction is a "d10" configuration, not a "d9s1" one. This is because Au atom prefers d9s1 configuration when isolated, but prefers d10 configuration in usual chemical environments, which has to be taken into consideration for a meaningful analysis.

4. Align

    Now, let's perform the orbital alignment between the real state (the real calculated electronic structure of the whole molecule) and the reference state (a "fake" state generated by combining all fragment). Make sure all the files in the following list are present in your working directory.

    **Input**
    - *Au4L4-nosymm.fchk*
    - *Au4L4-nosymm.47*
    - *Au4L4.id.rev*
    - *Au-d10.47*
    - *Au4L4-nosymm-L-1-p0.47*
    - *Au4L4-nosymm-L-2-p0.47*
    - *Au4L4-nosymm-L-3-p0.47*
    - *Au4L4-nosymm-L-4-p0.47*

    **Command**
    ```
    python FAMO.py Au4L4-nosymm.fchk AU4L4-NOSYMM.47 Au4L4.id.rev > Au4L4.famo
    ```

    **Output**
    - *Au4L4-nosymm_famo.fchk*
    This file stores a new set of occupied orbitals which as a whole is equivalent to the canonical MOs as directly computed in Gaussian. The new set of orbitals is formed according to a "maximal overlap" principle, so those with large eigenvalues are "maximally overlapped" with corresponding orbitals in the reference state, while those with small eigenvalues are "minimally overlapped". Note that these orbitals are all fully occupied, which makes this method different from density partitioning methods (such as NBO, AdNDP, etc.).

    - *Au4L4.famo*
    This file collects the output message of the alignment script. Most contents should be self-explanatory. Some important information are given below with more descriptions.

    ```
    NO. of electron pairs in the complex: ...
    NO. of electron pairs in fragments: ...
    ```
    These lines print the NO. of electron pairs (in closed-shell calculations, or NO. of electrons for each spin in open-shell calculations) in each fragment as well as the whole molecule. Note that the electron count has not to be conserved during the fragmentation.

    ```
    Inactive    space: ...
    Deformed    space: ...
    Transferred space: ...
    Unique      space: ...
    ```
    These liens give the NO. of orbitals in each space.
    If electrons are not conserved in fragmentation, for example, there are N0 occupied orbitals in the molecule but a total number of Nf occupied orbitals in the fragments, then there would be |N0-Nf| orbitals in the Unique space. For those orbitals that are both occupied in molecule and fragments, they can be almost unchanged (Inactive), somewhat deformed (Deformed), or completely different (Transferred). They will be collected into corresponding space by their overlaps (singular values). This categorization is done, however, by a pre-defined threshold and is merely for ease of analysis.

    ```
    T+U space diagonalized by Fock matrix.
    I+D space aligned to FMOs.
    ```
    These two lines tell the users how the orbitals are canonicalized. By default, the Transferred (T) space and Unique (U) space are diagonalized by Fock matrix, while the Inactive (I) space and Deformed (D) space aligned back to FMOs.This, however, is not the only usage. In the future version of FAMO code, options will be open for users to select.

    ```
    FAMO saved to file ..._famo.fchk
    ```
    For visualization of FAMOs, open this fchk with your favourite visualizer and plot relevant orbitals.

    ```
    Mulliken population analysis of FAMO in the basis of AO
    Fragment contribution to each FAMO:
                  Eigenvalue     1 Au        2 Au        3 Au       ...
         1    F1    0.997903    0.999270    0.001081    0.001081    ...
         2    F1    0.998251    1.003619   -0.000713   -0.000713    ...
         .     .           .           .           .           .
         .     .           .           .           .           .
         .     .           .           .           .           .
        72    F8    0.947447   -0.006695   -0.006695   -0.006695    ...
        73     U   -0.582759    0.254871    0.254871    0.254871    ...
      U_total                   0.254871    0.254871    0.254871    ...
       Total      145.999999   18.879730   18.879730   18.879730    ...
      Nuclear     148.000000   19.000000   19.000000   19.000000    ...
       Charge       2.000001    0.120270    0.120270    0.120270    ...
    ```
    This section performs a Mulliken population analysis of each FAMO and output how much each fragment contributes to each FAMO.
    The first column is the orbital id.
    The second column is a short identifier for readers to quickly identify the nature of orbitals. Depending on how orbitals are canonicalized in the last step, this identifier will vary along with the eigenvalues given in the 3rd column. In the default scenario, 'F1' denotes that this is a FAMO that is aligned back to FMO of fragment 1, and the following eigenvalue is the corresponding overlap. 'T/U' denotes that this a FAMO from Transferred/Unique space, and the following eigenvalue is the correpsonding energy since the T/U space is canonicalized by energy.
    The remaining columns are the contributions of each fragment to each FAMO using a Mulliken population analysis algorithm.
    The last a few lines in this section prints a summary of the Mulliken population analysis. "U_total" denotes the total population of orbitals in the Unique space, while "Total" indicates all orbitals. "Nuclear" denotes the nuclear charge of each fragment. Hence "Total" minus "Nuclear" are just the Mulliken charge of the fragments, as given in the "Charge" line.

    For spin-polarized systems, everything will be done twice, for alpha and beta space respectively.

## Relevant Publications
For more details of the method, see: Zhang J.-X., Sheong F. K., et al. to be submitted
For application in gold clusters, see: https://pubs.acs.org/doi/10.1021/acs.inorgchem.0c00649
For application in metallaaromatics, see: Hua Y., Zhang J.-X., et al. to be submitted


