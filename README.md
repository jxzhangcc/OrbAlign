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

In the following tutorial, I will use the cluster compound [Au4(PH3)4]2+ as an example. All mentioned files can be found in OrbAlign/test/

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
    - *AU4L4-nosymm.47*

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
    python autoSVDO.py Au4L4-nosymm.gjf Au4L4.id
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
    - *Au4L4-nosymm-L-1-p0.gjf*
    - *Au4L4-nosymm-L-2-p0.gjf*
    - *Au4L4-nosymm-L-3-p0.gjf*
    - *Au4L4-nosymm-L-4-p0.gjf*

    **Command**
    ```
    for f in Au4L4-nosymm-*gjf;  do g16 $f; done
    ```

    **Output**
    - *AU4L4-NOSYMM-L-1-P0.47*
    - *AU4L4-NOSYMM-L-2-P0.47*
    - *AU4L4-NOSYMM-L-3-P0.47*
    - *AU4L4-NOSYMM-L-4-P0.47*

    Among all output files from Gaussian, only the .47 files are useful in the following step. If you prefer your own .47 filenames when performing Gaussian calculations, change them back to these default names or make sure they are consistent with the filenames in the `Au4L4.id.rev` file so that the program can still properly recognize them.

4. Align

    Now, let's perform the orbital alignment between the real state (the real calculated electronic structure of the whole molecule) and the reference state (a "fake" state generated by combining all fragment). Make sure all the files in the following list are present in your working directory.

    **Input**
    - *Au4L4-nosymm.fchk*
    - *AU4L4-NOSYMM.47*
    - *Au4L4.id.rev*
    - *AU.47*
    - *AU4L4-NOSYMM-L-1-P0.47*
    - *AU4L4-NOSYMM-L-2-P0.47*
    - *AU4L4-NOSYMM-L-3-P0.47*
    - *AU4L4-NOSYMM-L-4-P0.47*

    **Command**
    ```
    python SVDO.py Au4L4-nosymm.fchk AU4L4-NOSYMM.47 Au4L4.id.rev > svdo.log
    ```

    **Output**
    - *Au4L4-nosymm_svdmo.fchk*
    This file stores a new set of occupied orbitals which as a whole is equivalent to the canonical MOs as directly computed in Gaussian. The new set of orbitals is formed according to a "maximal overlap" principle, so those with large eigenvalues are "maximally overlapped" with corresponding orbitals in the reference state, while those with small eigenvalues are "minimally overlapped". Note that these orbitals are all fully occupied, which makes this method different from density partitioning methods (such as NBO, AdNDP, etc.).

    - *Au4L4-nosymm_svdfo.fchk*
    The corresponding orbitals that are occupied in the reference state and are biorthogonal with the above orbitals given in the `Au4L4-nosymm_svdmo.fchk` file.

    - *svdo.log*
    This file collects the output message of the alignment script. Most contents should be self-explanatory. Some important information are given below with more descriptions.

    ```
    #electrons in each fragment:
    [9 9 9 9 9 9 9 9]
    #electrons in whole molecule:
    73
    ```
    These lines print the NO. of electron pairs (in closed-shell calculations, or NO. of electrons for each spin in open-shell calculations) in each fragment as well as the whole molecule. Note that the electron count is not preserved during the fragmentation. The whole cluster [Au4(PH3)4]2+ has 73 electron pairs, while the NO. of electron pairs of fragments sum up to 72. Hence in this case, there will one "unmatched" orbital with eigenvalue exactly equal to 0.

    ```
    Overlap:
    [1.    1.    1.    1.    1.    1.    1.    1.    1.    1.    1.    1.
     ...
     0.995 0.995 0.992 0.992 0.992 0.991 0.991 0.991 0.977 0.915 0.915 0.915]
    0 orbitals have singular values smaller than 1e-4 and are moved to the other group.
    ```
    This is the list of eigenvalues of the cross-overlap matrix between real state and reference state. An eigenvalue of 1 means that this orbital is completely identical in both states. An eigenvalue close to 1 means a small deformation and an eigenvalue close to 0 means significant deformation. In this example, most electrons remain unchanged because the electrons of [Au]+ and [PH3] ligands are either core electrons or inert bonding electrons. The four electron pairs involved in [PH3]->[Au]+ donation are the most active electrons among all, whose associated eigenvalues are still larger than 0.9, indicating that they are almost unchanged during the aggregation of all eight fragments. If some orbitals have eigenvalues smaller than a certain threshold (by default 1e-4), these orbitals will also be considered as "unmatched".

    ```
    Real orbitals:
    #orbitals to diagonalize by Fock matrix: (1, 380)
    Diagonalized Fock matrix:
    [-0.583]
    ```
    If there are orbitals appearing as "unmatched" during the alignment, then the OrbAlign program will further perform a diagonalization of Fock matrix within the space formed by the unmatch orbitals. The diagonalized Fock matrix can then be taken as "pseudo orbital energies" of the new orbitals.

    ```
    Reference orbitals:
    #orbitals to diagonalize by Fock matrix: (0, 380)
    Diagonalized Fock matrix
    []
    ```
    A similar diagonalization step is performed on the "unmatched" orbitals of the reference state if a) there are more electrons in fragments than in whole molecule; or b) there are very small eigenvalues appearing during alignment.

## Relevant Publications
For application in gold clusters, see: https://pubs.acs.org/doi/10.1021/acs.inorgchem.0c00649
For application in metallaaromatics, see: Hua Y., Zhang J.-X., et al, to be submitted


