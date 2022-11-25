Running the OLED DFT workflow
#############################

This is a minimal guide to running my DFT workflow. The first step is to install Pylada (see the Installation section of the other tutorial).

Setting up the calculations
===========================

The first thing you will need is a folder containing all the molecules in the POSCAR form. I convert the files using `Open Babel
<http://openbabel.org/wiki/Category:Installation>`_. but the resulting POSCARS have 0 unit cells so I use a different script to add a unit cell the the files. If you have obabel and Pylada installed, you can use the script `tovasp.sh` to do this automatically:

 .. code-block:: console

    bash tovasp.sh dir_containing_mols/

This will create a folder named ``POSCARs`` with the converted POSCARs. You need to make sure the file ``center.py`` is in the current directory.

On the supercomputer create a directory in ``scratch`` with the following files and folders:

 .. code-block:: console

    POSCARs #the directory we just created
    interp.py
    reposition.py
    stretch.py
    extract_data_b3lyp.py
    HT_workflow.ipy

In your ``home`` directory create a directory named ``pylada_chains`` add ``OLED_workflow_dryrun.py`` to it.

In the python environment where you have Pylada. You will have to install the following python packages:

 .. code-block:: console

    pip install ase pymatgen

Launching the calculations
==========================

To launch calculations go to the directory you just created in ``scratch`` and to this:

 .. code-block:: console

    ipython
    In [1]: %run HT_workflow.ipy
    In [2]: %launch scattered --account=rrg-ovoznyy --walltime=24:00:00 --ppn=40 --queue=compute

You can adjust the launch command according to the computer you are working on and the allocation you want to use. This will launch calculations for all the molecules in the ``POSCARs`` directory. If calculations hit the walltime, you can safely repeat this step, jobs that are finsihed will finsih as soon as they strat and won't be overwritten and jobs that did not fisnish will continue where they left off.

Extracting the results
======================

Once all calculations are finished use the following to get a text file containg all the results:

 .. code-block:: console

    python extract_data.py

This will create a text file named ``results.txt`` with the following structure:

+----------------+----------------+--------------+-------------------------------------+-----------------------------------+------------------------------+-----------------------+
| file name      | Energy Gap     | k\ :sub:`r`\ | Irridium Spin Orbit Coupling Energy |  Activation barrier (not relaxed) | Activation barrier (relaxed) | weakest bond strength | 
+----------------+----------------+--------------+-------------------------------------+-----------------------------------+------------------------------+-----------------------+

To obtain wavelength

.. math::

   \lambda = \frac{1240}{E_{gap}}

k\ :sub:`r`\  is consistantly overestimated, so I adjust it with

.. math::

   k_r = \frac{k_r}{15}

k\ :sub:`nr`\  is obtained with:

.. math::

   k_{nr} = k_{nr}(T) + k_{ISC} = Ae^{-\frac{E_{rel}}{B k_B T}} + Ce^{-E_{gap}} = 10^7e^{-\frac{E_{rel}}{3 \times 0.02585}} + 10^5e^{-E_{gap}}

   
Good luck!
