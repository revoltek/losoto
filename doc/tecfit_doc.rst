*****************************************
Obtaining Phase Solutions for TEC Fitting
*****************************************

TEC fitting requires phase solutions in multiple directions and bands. These solutions can be obtained in a number of ways using a number of calibration packages. Usage of the ``ion_peeling.py`` peeling script is described here.

The ``ion_peeling.py`` script can be used to obtain phase solutions towards multiple directions that are defined either as single sources or as patches of multiple sources. The script is intended to be run on a directory containing all measurement sets for all fields and bands. The script then uses the ``gsm.py`` script to identify a set of directions to use for peeling (a custom sky model may be specified instead). The script selects sources for peeling by estimating their apparent fluxes in each field and band and then determines which sources need to be peeled. Peeling is then performed as needed on all measurement sets and the solutions are collected in H5parm files for use with LoSoTo in the TECFIT operation.

The following options are available in the ``ionfit.py`` peeling script:

::

    Usage: ion_peeling.py <RA (deg)> <Dec (deg)> <radius (deg)> <outfile>

    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -i INDIR, --indir=INDIR
                            Input directory [default: .]
      -o OUTDIR, --outdir=OUTDIR
                            Output directory [default: Peeled]
      -f FLUXCUT, --fluxcut=FLUXCUT
                            Minimum apparent flux at 60 MHz in Jy for calibrators
                            [default: 15.0]
      -x FLUXBIN, --fluxbin=FLUXBIN
                            Target flux per bin at 60 MHz in Jy for tesselation
                            [default: 10.0]
      -g GSM, --gsm=GSM     Global sky model BBS file, to use instead of gsm.py
                            [default: none]
      -m MAJCUT, --majcut=MAJCUT
                            Maximum major axis size in arcmin for calibrators
                            [default: none]
      -b NBANDS, --nbands=NBANDS
                            Minimum number of bands that a calibrator must have to
                            be used [default: 8]
      -n NCORES, --ncores=NCORES
                            Number of simultaneous bands to calibrate [default: 8]
      -v, --verbose         Set verbose output and interactive mode [default:
                            False]
      -p, --patches         Use patches instead of single sources for calibrators?
                            [default: False]
      -t, --timecorr        Use time-correlated solutions? [default: False]
      -d, --dryrun          Do a dry run (no calibration is done) [default: False]
      -s SOLINT, --solint=SOLINT
                            Solution interval to use for phase solutions (# of
                            time slots) [default: 5]
      -c, --clobber         Clobber existing output files? [default: False]



******************
Fitting TEC Values
******************

The TECFIT operation fits phase solutions in multiple direction and bands to derive TEC values per station for each direction and time. The TECFIT operation uses the ``lofar.expion.baselinefitting.fit()`` function to fit a TEC value to the phases. The term that is minimized includes all baselines, so there is no preferred reference station, and the residual is computed as the complex number :math:`1.0 - exp(1i phasedifference)`, which is zero when the phase difference is a multiple of :math:`2\pi`.

Collecting Phase Solutions
--------------------------
The TECFIT operation scans all solution sets of the input H5parm file and searches for direction-dependent and direction-independent phase solutions. All such phase solutions are sorted by band (defined as all frequencies within the given tolerance that defaults to 1 MHz), direction, and field. Direction-independent phases will be added to the direction-dependent ones to remove field-to-field differences resulting from different direction-independent solutions. Therefore, the input H5parm file should contain only solutions to be used in TEC fitting. Solutions for each band and field should be saved to separate solution sets, as should direction-dependent and direction-independent solutions.

TECFIT Options
-------------------

The following options are available for the TECFIT operation:

.. glossary::

    Algorithm
        This option sets the algorithm used to obtain purely non-instrumental effects from the input phases. The ``sourcediff`` algorithm uses the difference of the phase solutions between sources to remove unknown instrumental effects. Currently, only the ``sourcediff`` algorithm is available.

    MinBands
        This option sets the minimum number of bands that a source must have to be used during fitting. Sources with phase solutions in fewer than ``MinBands`` number of bands are ignored during fitting.

    MaxStations
        This option sets the maximum number of stations that are used in the fitting. The ``MaxStations`` number of stations closest to the core will be used.

    OutSoltab
        This parameter sets the output solution table (of type ``tec``) in which the TEC values will be stored. This table is needed by the TECSCREEN operation when fitting a TEC screen.



*************************************
Fitting TEC Screens to the TEC values
*************************************

The TECSCREEN operation uses the TEC values derived by the TECFIT operation to fit TEC screens that can be used by the AWimager to correct for direction-dependent effects at any location in the image. Screens are derived assuming a single-layer ionosphere through a decomposition of the TEC values into KL basis vectors. Details are given in Intema et al. (2009).

TECSCREEN Options
-----------------

The following options are available for the TECSCREEN operation:

.. glossary::

    Height
        This option sets the height of the single-layer screen in meters.

    Order
        This options sets the maximum order of the basis vectors used in the screen. Higher orders allow the screen to trace finer structures, but such structures may not be reliable unless they are properly constrained by the pierce-point coverage.

    OutSoltab
        This parameter sets the output solution table (of type ``tecscreen``) in which the TEC values that define the screen are stored. This table is needed by ``H5parm_exporter.py`` when exporting the TEC screen values to a parmdb for use with the AWimager.


*********************************
Exporting TEC Screens to a ParmDB
*********************************

Once a TEC screen has been derived, it can be exported to a parmdb for use with BBS and the AWimager using the ``H5parm_exporter.py`` tool. The solution set containing the TEC values of the screen (in a ``tecscreen`` solution table) must be specified in the ``H5parm_exporter.py`` call.

If the ``sourcediff`` algorithm was used during TEC fitting, a direction-independent calibration should be done with BBS before imaging. This calibration is necessary to recover the overall normalization of the screen, which was lost during the TECFIT operation, and any remaining instrumental effects. The calibration should use the parmdb containing the TEC screen values and the BBS parset should enable the following options:

::

    Model.Ionosphere.Enable = T
    Model.Ionosphere.Type = EXPION

To use the TEC screen during imaging with the AWimager, the following options should be set in the ``awimager`` call:

::

    applyIonosphere = 1
    parmdbname = parmdb



