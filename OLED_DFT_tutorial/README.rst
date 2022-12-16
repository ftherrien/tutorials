Running the OLED DFT Workflow
#############################

This is a minimal guide to running my DFT workflow. The first step is to install Pylada (see the Installation section of the other tutorial).

Setting up the calculations
===========================

Making the input files
----------------------

The first thing you will need is a folder containing all the molecules in the POSCAR form. I convert the files using `Open Babel
<http://openbabel.org/wiki/Category:Installation>`_. but the resulting POSCARS have 0 unit cells so I use a different script to add a unit cell to the files. If you have obabel and Pylada installed, you can use the script `tovasp.sh` to do this automatically:

.. code-block::
    
   bash tovasp.sh dir_containing_mols/

This will create a folder named ``POSCARs`` with the converted POSCARs. You need to make sure the file ``center.py`` is in the current directory.

.. note::

   If you want the scripts to work out of the box, the POSCARs need to be named in the form ``POSCAR_[entry]``.


Run directory
-------------

On the supercomputer create a directory in ``scratch`` with the following files and folders:

.. code-block::
    
   POSCARs #the directory we just created
   mols #a directory containing the equivalent entries in mol format
   interp.py
   reposition.py
   stretch.py
   extract_images.py
   extract_data.py
   mcu.py
   HT_workflow.ipy

Open ``HT_workflow.ipy``, at lines 29-30 set the location of the VASP executable and the pseudo-potentials:

.. code-block:: python
    
    ############### setting up the functional
    vasp=Relax()
    vasp.has_nlep = False
    
    vasp.program = '/scinet/niagara/software/commercial/vasp-6.1.2/bin/vasp_ncl' #<- this line
    pseudoDir = '/scinet/niagara/software/commercial/vasp-6.1.2/potpaw_PBE/' #<- and this line

    vasp.add_specie = "C", pseudoDir + "/C"
    vasp.add_specie = "H", pseudoDir + "/H"

.. note::

   You need to specify the path to ``vasp_ncl`` because spin orbit coupling requires non-collinear calculations.

Then, in your ``home`` directory create a directory named ``pylada_chains`` and add ``OLED_workflow_dryrun.py`` to it. Then add the following to your ``~/.bashrc``:

.. code-block::
   
   export PYTHONPATH=:$PYTHONPATH:$HOME/pylada_chains
   export VASP_PP_PATH="/scinet/niagara/software/commercial/vasp-6.1.2/" # <-- Path to vasp for ASE

The second line contains the path to VASP for ASE to read. Don't forget to source your ``~/.bashrc``:

.. code-block::
   
   source ~/.bashrc

Additional required packages
----------------------------

In the python environment where you have installed Pylada. You will have to install the following python packages:

.. code-block::
    
   pip install ase pymatgen rdkit

    
Launching the calculations
==========================

To launch calculations go to the directory you just created in ``scratch`` and to this:

.. code-block::
    
   ipython
   In [1]: %run HT_workflow.ipy
   In [2]: %launch scattered --account=rrg-ovoznyy --walltime=24:00:00 --ppn=40 --queue=compute

You can adjust the launch command according to the computer you are working on and the allocation you want to use. This will launch calculations for all the molecules in the ``POSCARs`` directory. If calculations hit the wall time, you can safely repeat this step, jobs that are finished will finish as soon as they start and won't be overwritten and jobs that did not finish will continue where they left off.

Extracting the results
======================

Once all calculations are finished, use the following to get a text file containing all the results:

.. code-block::
    
   python extract_data.py

This will create a text file named ``results.txt`` with the following structure:

+-----------+------------+-------------------------+----------------------+----------------------------+-----------+-------------+
| file name | Wavelength | PLQY unrelaxed barriers | PLQY relaxed barrier | PLQY fully relaxed barrier | Weak bond |  Not N bond |
+-----------+------------+-------------------------+----------------------+----------------------------+-----------+-------------+
|           |            | Only if blue            | Only if blue         | Only if blue               | if weak   |  if not N   |
+-----------+------------+-------------------------+----------------------+----------------------------+-----------+-------------+


.. math::

   \lambda = \frac{1240}{E_{gap}}

k\ :sub:`r`\  is consistently overestimated, so I adjust it with

.. math::

   k_r = \frac{k_r}{15}

k\ :sub:`nr`\  is obtained with:

.. math::

   k_{nr} = k_{nr}(T) + k_{ISC} = Ae^{-\frac{E_{rel}}{B k_B T}} + Ce^{-E_{gap}} = 10^7e^{-\frac{E_{rel}}{3 \times 0.02585}} + 10^5e^{-E_{gap}}


"PLQY unrelaxed barriers" uses barriers found by manually stretching the bonds by sets amounts and stopping when the LUMO has high metal character. This is the original method.

"PLQY relaxed barrier" takes the structure with the lowest barrier from the previous step and relaxes it while keeping all Ir bonds fixed.

"PLQY relaxed barrier" starts from the previous step and keeps only the stertched bond fixed while relaxing all other bonds.
   
Good luck!
