#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation implements screen plotting

from losoto.lib_operations import *
from losoto.operations.directionscreen import _calc_piercepoint
from losoto._logging import logger as logging

logging.debug('Loading PLOTSCREEN module.')

def _run_parser(soltab, parser, step):
    resSoltab = parser.getstr( step, "resSoltab", '' )
    minZ, maxZ = parser.getarrayfloat( step, "MinMax", [0.0, 0.0] )
    prefix = parser.getstr( step, "Prefix", '' )
    remove_gradient = parser.getbool( step, "RemoveGradient", False )
    show_source_names = parser.getbool( step, "ShowSourceNames", False )
    ncpu = parser.getint( step, "npcu", 0 )

    parser.checkSpelling( step, soltab, ['resSoltab', 'minZ', 'maxZ', 'prefix', 'remove_gradient', 'show_source_names'])
    return run(soltab, resSoltab, minZ, maxZ, prefix, remove_gradient, show_source_names, ncpu)


def _phase_cm():
    """
    Returns colormap for phase plots
    """
    from matplotlib.colors import ListedColormap

    cm_data = [[ 0.65830839, 0.46993917, 0.04941288],
               [ 0.66433742, 0.4662019 , 0.05766473],
               [ 0.67020869, 0.46248014, 0.0653456 ],
               [ 0.67604299, 0.45869838, 0.07273174],
               [ 0.68175228, 0.45491407, 0.07979262],
               [ 0.6874028 , 0.45108417, 0.08667103],
               [ 0.6929505 , 0.44723893, 0.09335869],
               [ 0.69842619, 0.44335768, 0.09992839],
               [ 0.7038123 , 0.43945328, 0.1063871 ],
               [ 0.70912069, 0.43551765, 0.11277174],
               [ 0.71434524, 0.43155576, 0.11909348],
               [ 0.71949289, 0.42756272, 0.12537606],
               [ 0.72455619, 0.4235447 , 0.13162325],
               [ 0.72954895, 0.41949098, 0.13786305],
               [ 0.73445172, 0.41541774, 0.14408039],
               [ 0.73929496, 0.41129973, 0.15032217],
               [ 0.74403834, 0.40717158, 0.15654335],
               [ 0.74873695, 0.40298519, 0.16282282],
               [ 0.75332319, 0.39880107, 0.16907566],
               [ 0.75788083, 0.39454245, 0.17542179],
               [ 0.7623326 , 0.39028096, 0.18175915],
               [ 0.76673205, 0.38596549, 0.18816819],
               [ 0.77105247, 0.38162141, 0.19461532],
               [ 0.77529528, 0.37724732, 0.20110652],
               [ 0.77948666, 0.37281509, 0.2076873 ],
               [ 0.78358534, 0.36836772, 0.21429736],
               [ 0.78763763, 0.363854  , 0.22101648],
               [ 0.79161134, 0.35930804, 0.2277974 ],
               [ 0.79550606, 0.3547299 , 0.23464353],
               [ 0.79935398, 0.35007959, 0.24161832],
               [ 0.80311671, 0.34540152, 0.24865892],
               [ 0.80681033, 0.34067452, 0.25580075],
               [ 0.8104452 , 0.33588248, 0.26307222],
               [ 0.8139968 , 0.33105538, 0.27043183],
               [ 0.81747689, 0.32617526, 0.27791096],
               [ 0.82089415, 0.32122629, 0.28553846],
               [ 0.82422713, 0.3162362 , 0.29327617],
               [ 0.82747661, 0.31120154, 0.30113388],
               [ 0.83066399, 0.30608459, 0.30917579],
               [ 0.83376307, 0.30092244, 0.31734921],
               [ 0.83677286, 0.29571346, 0.32566199],
               [ 0.83969693, 0.29044723, 0.33413665],
               [ 0.84253873, 0.28511151, 0.34279962],
               [ 0.84528297, 0.27972917, 0.35162078],
               [ 0.84792704, 0.27430045, 0.36060681],
               [ 0.85046793, 0.26882624, 0.36976395],
               [ 0.85291056, 0.26328859, 0.37913116],
               [ 0.855242  , 0.25770888, 0.38868217],
               [ 0.85745673, 0.25209367, 0.39841601],
               [ 0.85955023, 0.24644737, 0.40833625],
               [ 0.86151767, 0.24077563, 0.41844557],
               [ 0.86335392, 0.23508521, 0.42874606],
               [ 0.86505685, 0.22937288, 0.43926008],
               [ 0.86661606, 0.22366308, 0.44996127],
               [ 0.86802578, 0.21796785, 0.46084758],
               [ 0.86928003, 0.21230132, 0.47191554],
               [ 0.87037274, 0.20667988, 0.48316015],
               [ 0.87129781, 0.2011224 , 0.49457479],
               [ 0.87204914, 0.19565041, 0.50615118],
               [ 0.87262076, 0.19028829, 0.51787932],
               [ 0.87300686, 0.18506334, 0.5297475 ],
               [ 0.8732019 , 0.18000588, 0.54174232],
               [ 0.87320066, 0.1751492 , 0.55384874],
               [ 0.87299833, 0.17052942, 0.56605016],
               [ 0.87259058, 0.16618514, 0.57832856],
               [ 0.87197361, 0.16215698, 0.59066466],
               [ 0.87114414, 0.15848667, 0.60303881],
               [ 0.87009966, 0.15521687, 0.61542844],
               [ 0.86883823, 0.15238892, 0.62781175],
               [ 0.86735858, 0.15004199, 0.64016651],
               [ 0.8656601 , 0.14821149, 0.65247022],
               [ 0.86374282, 0.14692762, 0.66470043],
               [ 0.86160744, 0.14621386, 0.67683495],
               [ 0.85925523, 0.14608582, 0.68885204],
               [ 0.85668805, 0.14655046, 0.70073065],
               [ 0.85390829, 0.14760576, 0.71245054],
               [ 0.85091881, 0.14924094, 0.7239925 ],
               [ 0.84772287, 0.15143717, 0.73533849],
               [ 0.84432409, 0.15416865, 0.74647174],
               [ 0.84072639, 0.15740403, 0.75737678],
               [ 0.83693394, 0.16110786, 0.76803952],
               [ 0.83295108, 0.16524205, 0.77844723],
               [ 0.82878232, 0.16976729, 0.78858858],
               [ 0.82443225, 0.17464414, 0.7984536 ],
               [ 0.81990551, 0.179834  , 0.80803365],
               [ 0.81520674, 0.18529984, 0.8173214 ],
               [ 0.81034059, 0.19100664, 0.82631073],
               [ 0.80531176, 0.1969216 , 0.83499645],
               [ 0.80012467, 0.20301465, 0.84337486],
               [ 0.79478367, 0.20925826, 0.8514432 ],
               [ 0.78929302, 0.21562737, 0.85919957],
               [ 0.78365681, 0.22209936, 0.86664294],
               [ 0.77787898, 0.22865386, 0.87377308],
               [ 0.7719633 , 0.23527265, 0.88059043],
               [ 0.76591335, 0.24193947, 0.88709606],
               [ 0.7597325 , 0.24863985, 0.89329158],
               [ 0.75342394, 0.25536094, 0.89917908],
               [ 0.74699063, 0.26209137, 0.90476105],
               [ 0.74043533, 0.2688211 , 0.91004033],
               [ 0.73376055, 0.27554128, 0.91502   ],
               [ 0.72696862, 0.28224415, 0.91970339],
               [ 0.7200616 , 0.2889229 , 0.92409395],
               [ 0.71304134, 0.29557159, 0.92819525],
               [ 0.70590945, 0.30218508, 0.9320109 ],
               [ 0.69866732, 0.30875887, 0.93554451],
               [ 0.69131609, 0.31528914, 0.93879964],
               [ 0.68385669, 0.32177259, 0.94177976],
               [ 0.6762898 , 0.32820641, 0.94448822],
               [ 0.6686159 , 0.33458824, 0.94692818],
               [ 0.66083524, 0.3409161 , 0.94910264],
               [ 0.65294785, 0.34718834, 0.95101432],
               [ 0.64495358, 0.35340362, 0.95266571],
               [ 0.63685208, 0.35956083, 0.954059  ],
               [ 0.62864284, 0.3656591 , 0.95519608],
               [ 0.62032517, 0.3716977 , 0.95607853],
               [ 0.61189825, 0.37767607, 0.95670757],
               [ 0.60336117, 0.38359374, 0.95708408],
               [ 0.59471291, 0.3894503 , 0.95720861],
               [ 0.58595242, 0.39524541, 0.95708134],
               [ 0.5770786 , 0.40097871, 0.95670212],
               [ 0.56809041, 0.40664983, 0.95607045],
               [ 0.55898686, 0.41225834, 0.95518556],
               [ 0.54976709, 0.41780374, 0.95404636],
               [ 0.5404304 , 0.42328541, 0.95265153],
               [ 0.53097635, 0.42870263, 0.95099953],
               [ 0.52140479, 0.43405447, 0.94908866],
               [ 0.51171597, 0.43933988, 0.94691713],
               [ 0.50191056, 0.44455757, 0.94448311],
               [ 0.49198981, 0.44970607, 0.94178481],
               [ 0.48195555, 0.45478367, 0.93882055],
               [ 0.47181035, 0.45978843, 0.93558888],
               [ 0.46155756, 0.46471821, 0.93208866],
               [ 0.45119801, 0.46957218, 0.92831786],
               [ 0.44073852, 0.47434688, 0.92427669],
               [ 0.43018722, 0.47903864, 0.9199662 ],
               [ 0.41955166, 0.4836444 , 0.91538759],
               [ 0.40884063, 0.48816094, 0.91054293],
               [ 0.39806421, 0.49258494, 0.90543523],
               [ 0.38723377, 0.49691301, 0.90006852],
               [ 0.37636206, 0.50114173, 0.89444794],
               [ 0.36546127, 0.5052684 , 0.88857877],
               [ 0.35454654, 0.5092898 , 0.88246819],
               [ 0.34363779, 0.51320158, 0.87612664],
               [ 0.33275309, 0.51700082, 0.86956409],
               [ 0.32191166, 0.52068487, 0.86279166],
               [ 0.31113372, 0.52425144, 0.85582152],
               [ 0.3004404 , 0.52769862, 0.84866679],
               [ 0.28985326, 0.53102505, 0.84134123],
               [ 0.27939616, 0.53422931, 0.83386051],
               [ 0.26909181, 0.53731099, 0.82623984],
               [ 0.258963  , 0.5402702 , 0.81849475],
               [ 0.24903239, 0.54310763, 0.8106409 ],
               [ 0.23932229, 0.54582448, 0.80269392],
               [ 0.22985664, 0.54842189, 0.79467122],
               [ 0.2206551 , 0.55090241, 0.78658706],
               [ 0.21173641, 0.55326901, 0.77845533],
               [ 0.20311843, 0.55552489, 0.77028973],
               [ 0.1948172 , 0.55767365, 0.76210318],
               [ 0.1868466 , 0.55971922, 0.75390763],
               [ 0.17921799, 0.56166586, 0.74571407],
               [ 0.1719422 , 0.56351747, 0.73753498],
               [ 0.16502295, 0.56527915, 0.72937754],
               [ 0.15846116, 0.566956  , 0.72124819],
               [ 0.15225499, 0.56855297, 0.71315321],
               [ 0.14639876, 0.57007506, 0.70509769],
               [ 0.14088284, 0.57152729, 0.69708554],
               [ 0.13569366, 0.57291467, 0.68911948],
               [ 0.13081385, 0.57424211, 0.68120108],
               [ 0.12622247, 0.57551447, 0.67333078],
               [ 0.12189539, 0.57673644, 0.66550792],
               [ 0.11780654, 0.57791235, 0.65773233],
               [ 0.11392613, 0.5790468 , 0.64999984],
               [ 0.11022348, 0.58014398, 0.64230637],
               [ 0.10666732, 0.58120782, 0.63464733],
               [ 0.10322631, 0.58224198, 0.62701729],
               [ 0.0998697 , 0.58324982, 0.61941001],
               [ 0.09656813, 0.58423445, 0.61181853],
               [ 0.09329429, 0.58519864, 0.60423523],
               [ 0.09002364, 0.58614483, 0.5966519 ],
               [ 0.08673514, 0.58707512, 0.58905979],
               [ 0.08341199, 0.58799127, 0.58144971],
               [ 0.08004245, 0.58889466, 0.57381211],
               [ 0.07662083, 0.58978633, 0.56613714],
               [ 0.07314852, 0.59066692, 0.55841474],
               [ 0.06963541, 0.5915367 , 0.55063471],
               [ 0.06610144, 0.59239556, 0.54278681],
               [ 0.06257861, 0.59324304, 0.53486082],
               [ 0.05911304, 0.59407833, 0.52684614],
               [ 0.05576765, 0.5949003 , 0.5187322 ],
               [ 0.05262511, 0.59570732, 0.51050978],
               [ 0.04978881, 0.5964975 , 0.50216936],
               [ 0.04738319, 0.59726862, 0.49370174],
               [ 0.04555067, 0.59801813, 0.48509809],
               [ 0.04444396, 0.59874316, 0.47635   ],
               [ 0.04421323, 0.59944056, 0.46744951],
               [ 0.04498918, 0.60010687, 0.45838913],
               [ 0.04686604, 0.60073837, 0.44916187],
               [ 0.04988979, 0.60133103, 0.43976125],
               [ 0.05405573, 0.60188055, 0.4301812 ],
               [ 0.05932209, 0.60238289, 0.42040543],
               [ 0.06560774, 0.60283258, 0.41043772],
               [ 0.07281962, 0.60322442, 0.40027363],
               [ 0.08086177, 0.60355283, 0.38990941],
               [ 0.08964366, 0.60381194, 0.37934208],
               [ 0.09908952, 0.60399554, 0.36856412],
               [ 0.10914617, 0.60409695, 0.35755799],
               [ 0.11974119, 0.60410858, 0.34634096],
               [ 0.13082746, 0.6040228 , 0.33491416],
               [ 0.14238003, 0.60383119, 0.323267  ],
               [ 0.1543847 , 0.60352425, 0.31138823],
               [ 0.16679093, 0.60309301, 0.29931029],
               [ 0.17959757, 0.60252668, 0.2870237 ],
               [ 0.19279966, 0.60181364, 0.27452964],
               [ 0.20634465, 0.60094466, 0.2618794 ],
               [ 0.22027287, 0.5999043 , 0.24904251],
               [ 0.23449833, 0.59868591, 0.23611022],
               [ 0.24904416, 0.5972746 , 0.2230778 ],
               [ 0.26382006, 0.59566656, 0.21004673],
               [ 0.2788104 , 0.5938521 , 0.19705484],
               [ 0.29391494, 0.59183348, 0.18421621],
               [ 0.3090634 , 0.58961302, 0.17161942],
               [ 0.32415577, 0.58720132, 0.15937753],
               [ 0.3391059 , 0.58461164, 0.14759012],
               [ 0.35379624, 0.58186793, 0.13637734],
               [ 0.36817905, 0.5789861 , 0.12580054],
               [ 0.38215966, 0.57599512, 0.1159504 ],
               [ 0.39572824, 0.57290928, 0.10685038],
               [ 0.40881926, 0.56975727, 0.09855521],
               [ 0.42148106, 0.56654159, 0.09104002],
               [ 0.43364953, 0.56329296, 0.08434116],
               [ 0.44538908, 0.56000859, 0.07841305],
               [ 0.45672421, 0.5566943 , 0.07322913],
               [ 0.46765017, 0.55336373, 0.06876762],
               [ 0.47819138, 0.5500213 , 0.06498436],
               [ 0.48839686, 0.54666195, 0.06182163],
               [ 0.49828924, 0.5432874 , 0.05922726],
               [ 0.50789114, 0.53989827, 0.05714466],
               [ 0.51722475, 0.53649429, 0.05551476],
               [ 0.5263115 , 0.53307443, 0.05427793],
               [ 0.53517186, 0.52963707, 0.05337567],
               [ 0.54382515, 0.52618009, 0.05275208],
               [ 0.55228947, 0.52270103, 0.05235479],
               [ 0.56058163, 0.51919713, 0.0521356 ],
               [ 0.56871719, 0.51566545, 0.05205062],
               [ 0.57671045, 0.51210292, 0.0520602 ],
               [ 0.5845745 , 0.50850636, 0.05212851],
               [ 0.59232129, 0.50487256, 0.05222299],
               [ 0.5999617 , 0.50119827, 0.05231367],
               [ 0.60750568, 0.49748022, 0.05237234],
               [ 0.61496232, 0.49371512, 0.05237168],
               [ 0.62233999, 0.48989963, 0.05228423],
               [ 0.62964652, 0.48603032, 0.05208127],
               [ 0.63688935, 0.48210362, 0.05173155],
               [ 0.64407572, 0.4781157 , 0.0511996 ],
               [ 0.65121289, 0.47406244, 0.05044367],
               [ 0.65830839, 0.46993917, 0.04941288]]
    cm = ListedColormap(cm_data, name=__file__)

    return cm


def _calculate_screen(inscreen, residuals, pp, N_piercepoints, k, east, north, up,
    T, Nx, Ny, sindx, height, beta_val, r_0, is_phase, outQueue):
    """
    Calculates screen images

    Parameters
    ----------
    inscreen : array
        Array of screen values at the piercepoints
    residuals : array
        Array of screen residuals at the piercepoints
    pp : array
        Array of piercepoint locations
    N_piercepoints : int
        Number of pierce points
    k : int
        Time index
    east : array
        East array
    north : array
        North array
    up : array
        Up array
    T : array
        T array
    Nx : int
        Number of pixels in x for screen
    Ny : int
        Number of pixels in y for screen
    sindx : int
        Station index
    height : float
        height of screen (m)
    beta_val : float
        power-law index for phase structure function (5/3 =>
        pure Kolmogorov turbulence)
    r_0 : float
        scale size of phase fluctuations
    is_phase : bool
        input screen is a phase screen

    """
    from numpy import kron, concatenate, newaxis
    from numpy.linalg import pinv, norm
    import numpy as np

    screen = np.zeros((Nx, Ny))

    if height == 0.0:
        pp1 = pp[:, :]
    else:
        pp1 = np.dot(pp[:, :], T)

    min_xy = np.amin(pp1, axis=0)
    max_xy = np.amax(pp1, axis=0)
    extent = max_xy - min_xy
    lowerk = min_xy - 0.1 * extent
    upperk = max_xy + 0.1 * extent
    im_extent_mk = upperk - lowerk
    pix_per_mk = Nx / im_extent_mk[0]
    m_per_pixk = 1.0 / pix_per_mk

    xr = np.arange(lowerk[0], upperk[0], m_per_pixk)
    yr = np.arange(lowerk[1], upperk[1], m_per_pixk)
    D = np.resize(pp, (N_piercepoints, N_piercepoints, 3))
    D = np.transpose(D, (1, 0, 2)) - D
    D2 = np.sum(D**2, axis=2)
    C = -(D2 / r_0**2)**(beta_val / 2.0) / 2.0
    f = inscreen.reshape(N_piercepoints)
    fitted_tec = np.dot(C, f) + residuals
    if is_phase:
        fitted_tec = normalize_phase(fitted_tec)
    for i, xi in enumerate(xr[0: Nx]):
        for j, yi in enumerate(yr[0: Ny]):
            if height == 0.0:
                p = np.array([xi, yi, 0.0])
            else:
                p, airmass = _calc_piercepoint(np.dot(np.array([xi, yi]), np.array([east, north])), up, height)
            d2 = np.sum(np.square(pp - p), axis=1)
            c = -(d2 / ( r_0**2 ))**(beta_val / 2.0) / 2.0
            screen[i, j] = np.dot(c, f)

    # Calculate the piercepoint coords, normalized to 1
    if height == 0.0:
        x = pp1[:, 0]
        y = pp1[:, 1]
    else:
        x = (pp1[:, 0] / 1000.0 - min_xy[0] / 1000.0) / extent[0] # km
        y = (pp1[:, 1] / 1000.0 - min_xy[1] / 1000.0) / extent[1] # km

    outQueue.put([k, fitted_tec, screen, x, y])


def _plot_frame(screen, fitted_phase1, residuals, weights, x, y, k, lower,
    upper, vmin, vmax, source_names, show_source_names, station_names, sindx,
    root_dir, prestr, is_image_plane,  midRA, midDec, order, is_phase, outQueue):
    """
    Plots screen images

    Parameters
    ----------
    screen : array
        Image of screen values
    fitted_phase1 : array
        Array of fitted phase values
    residuals : array
        Array of phase residuals at the piercepoints
    weights : array
        Array of weights at the piercepoints
    x : array
        Array of piercepoint x locations
    y : array
        Array of piercepoint y locations
    k : int
        Time index
    lower : array
        Array of lower limits for plot
    upper : array
        Array of upper limits for plot
    vmin : float
        minimum value for plot range
    vmax : float
        maximum value for plot range
    source_names : list
        List of source (direction) names
    show_source_names : bool
        label sources on screen plots
    order : int
        order of screen
    is_phase : bool
        True if screens are phase screens

    """
    if not 'matplotlib' in sys.modules:
        import matplotlib as mpl
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
        mpl.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt # after setting "Agg" to speed up
    from losoto.operations.stationscreen import _makeWCS, _circ_chi2
    import numpy as np
    try:
        try:
            from astropy.visualization.wcsaxes import WCSAxes
            hasWCSaxes = True
        except:
            from wcsaxes import WCSAxes
            hasWCSaxes = True
    except:
        hasWCSaxes = False
    from matplotlib.colors import LinearSegmentedColormap


    fig = plt.figure(figsize=(6,6))

    # Set colormap
    if is_phase:
        cmap = _phase_cm()
    else:
        cmap = plt.cm.jet
    sm = plt.cm.ScalarMappable(cmap=cmap,
        norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []

    if is_image_plane and hasWCSaxes:
        wcs = _makeWCS(midRA, midDec)
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
    else:
        plt.gca().set_aspect('equal')
        ax = plt.gca()

    s = []
    c = []
    xf = []
    yf = []
    weights = np.array(weights, dtype=float)
    nonflagged = np.where(weights > 0.0)
    for j in range(fitted_phase1.shape[0]):
        if weights[j] > 0.0:
            s.append(max(20, 200*np.sqrt(weights[j]/np.median(weights[nonflagged]))))
        else:
            s.append(120)
            xf.append(x[j])
            yf.append(y[j])
        c.append(sm.to_rgba(fitted_phase1[j]))

    if is_image_plane:
        min_x = np.min(x)
        max_x = np.max(x)
        min_y = np.min(y)
        max_y = np.max(y)
        extent_x = max_x - min_x
        extent_y = max_y - min_y
        lower = [min_x - 0.1 * extent_x, min_y - 0.1 * extent_y]
        upper = [max_x + 0.1 * extent_x, max_y + 0.1 * extent_y]
        Nx = screen.shape[0]
        pix_per_m = Nx / (upper[0] - lower[0])
        m_per_pix = 1.0 / pix_per_m
        xr = np.arange(lower[0], upper[0], m_per_pix)
        yr = np.arange(lower[1], upper[1], m_per_pix)
        lower = np.array([xr[0], yr[0]])
        upper = np.array([xr[-1], yr[-1]])
    else:
        # convert from m to km
        lower /= 1000.0
        upper /= 1000.0

    im = ax.imshow(screen.transpose([1, 0])[:, :],
        cmap = cmap,
        origin = 'lower',
        interpolation = 'nearest',
        extent = (lower[0], upper[0], lower[1], upper[1]),
        vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im)
    cbar.set_label('Value', rotation=270)

    ax.scatter(np.array(x), np.array(y), s=np.array(s), c=np.array(c), alpha=0.7, cmap=cmap, vmin=vmin, vmax=vmax, edgecolor='black')
    if len(xf) > 0:
        ax.scatter(xf, yf, s=120, c='k', marker='x')
    if show_source_names:
        labels = source_names
        for label, xl, yl in zip(labels, x, y):
            plt.annotate(
                label,
                xy = (xl, yl), xytext = (-2, 2),
                textcoords = 'offset points', ha = 'right', va = 'bottom')

    nsrcs = np.where(weights > 0.0)[0].size
    if is_phase:
        redchi2 =  _circ_chi2(residuals, weights) / (nsrcs-order)
    else:
        redchi2 =  np.sum(np.square(residuals) * weights) / (nsrcs-order)
    if sindx >= 0:
        plt.title('Station {0}, Time {1} (red. chi2 = {2:0.3f})'.format(station_names[sindx], k, redchi2))
    else:
        plt.title('Time {0}'.format(k))
    if is_image_plane:
        ax.set_xlim(lower[0], upper[0])
        ax.set_ylim(lower[1], upper[1])
        ax.set_aspect('equal')
        if hasWCSaxes:
            RAAxis = ax.coords['ra']
            RAAxis.set_axislabel('RA', minpad=0.75)
            RAAxis.set_major_formatter('hh:mm:ss')
            DecAxis = ax.coords['dec']
            DecAxis.set_axislabel('Dec', minpad=0.75)
            DecAxis.set_major_formatter('dd:mm:ss')
            ax.coords.grid(color='black', alpha=0.5, linestyle='solid')
            plt.xlabel("RA")
            plt.ylabel("Dec")
        else:
            plt.xlabel("RA (arb. units)")
            plt.ylabel("Dec (arb. units)")
    else:
        # Reverse the axis so that RA coord increases to left
        plt.xlim(upper[0], lower[0])
        plt.ylim(lower[1], upper[1])
        plt.xlabel('Projected Distance East-West (km)')
        plt.ylabel('Projected Distance North-South (km)')
    if sindx >= 0:
        plt.savefig(root_dir + '/' + prestr + '_station%0.4i' % sindx + '_frame%0.4i.png' % k, bbox_inches='tight')
    else:
        plt.savefig(root_dir + '/' + prestr + '_frame%0.4i.png' % k, bbox_inches='tight')
    plt.close(fig)


def _make_screen_plots(pp, inscreen, inresiduals, weights, station_names,
    station_positions, source_names, times, height, station_order, beta_val,
    r_0, prefix='frame_', remove_gradient=True, show_source_names=False,
    min_val=None, max_val=None, is_phase=False, midRA=0.0, midDec=0.0, ncpu=0):
    """
    Makes plots of screens

    Parameters
    ----------
    pp : array
        Array of piercepoint locations
    inscreen : array
        Array of screen values at the piercepoints with order [dir, time, ant]
    residuals : array
        Array of screen residuals at the piercepoints with order [dir, time, ant]
    weights : array
        Array of weights for each piercepoint with order [dir, time, ant]
    source_names: array
        Array of source names
    times : array
        Array of times
    height : float
        Height of screen (m)
    station_order : list
        List of order of screens per station (e.g., number of KL base vectors to keep)
    r_0 : float
        Scale size of phase fluctuations
    beta_val : float
        Power-law index for phase structure function
    prefix : str
        Prefix for output file names
    remove_gradient : bool
        Fit and remove a gradient from each screen
    show_source_names : bool
        Label sources on screen plots
    min_val : float
        Minimum value for plot range
    max_val : float
        Maximum value for plot range
    is_phase : bool
        Input screen is a phase screen
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees
    ncpu : int
        Number of CPUs to use

    """
    from numpy import kron, concatenate, newaxis
    from numpy.linalg import pinv, norm
    import numpy as np
    import os

    # avoids error if re-setting "agg" a second run of plot
    if not 'matplotlib' in sys.modules:
        import matplotlib as mpl
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
        mpl.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt # after setting "Agg" to speed up
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    # input check
    root_dir = os.path.dirname(prefix)
    if root_dir == '':
        root_dir = './'
    prestr = os.path.basename(prefix) + 'screen'
    try:
        os.makedirs(root_dir)
    except OSError:
        pass

    if height == 0.0:
        N_stations = 1 # screens are single-station screens
    else:
        N_stations = len(station_names)
    N_sources = len(source_names)
    N_times = len(times)
    N_piercepoints = N_sources * N_stations
    xp, yp, zp = station_positions[0, :] # use first station
    east = np.array([-yp, xp, 0])
    east = east / norm(east)

    north = np.array([-xp, -yp, (xp*xp + yp*yp)/zp])
    north = north / norm(north)

    up = np.array([xp, yp, zp])
    up = up / norm(up)

    T = concatenate([east[:, newaxis], north[:, newaxis]], axis=1)

    # Use pierce point locations of first and last time slots to estimate
    # required size of plot in meters
    if height == 0.0:
        is_image_plane = True # pierce points are image plane coords
        pp1_0 = pp[:, 0:2]
        pp1_1 = pp[:, 0:2]
    else:
        is_image_plane = False
        pp1_0 = np.dot(pp[0, :, :], T)
        pp1_1 = np.dot(pp[-1, :, :], T)

    max_xy = np.amax(pp1_0, axis=0) - np.amin(pp1_0, axis=0)
    max_xy_1 = np.amax(pp1_1, axis=0) - np.amin(pp1_1, axis=0)
    if max_xy_1[0] > max_xy[0]:
        max_xy[0] = max_xy_1[0]
    if max_xy_1[1] > max_xy[1]:
        max_xy[1] = max_xy_1[1]

    min_xy = np.array([0.0, 0.0])
    extent = max_xy - min_xy
    lower = min_xy - 0.1 * extent
    upper = max_xy + 0.1 * extent
    im_extent_m = upper - lower
    fitted_phase1 = np.zeros((N_piercepoints, N_times))

    Nx = 60 # set approximate number of pixels in screen
    pix_per_m = Nx / im_extent_m[0]
    m_per_pix = 1.0 / pix_per_m
    xr = np.arange(lower[0], upper[0], m_per_pix)
    yr = np.arange(lower[1], upper[1], m_per_pix)
    Nx = len(xr)
    Ny = len(yr)
    lower = np.array([xr[0], yr[0]])
    upper = np.array([xr[-1], yr[-1]])

    x = np.zeros((N_times, N_piercepoints)) # plot x pos of piercepoints
    y = np.zeros((N_times, N_piercepoints)) # plot y pos of piercepoints
    screen = np.zeros((Nx, Ny, N_times))

    if height == 0.0:
        for sindx in range(station_positions.shape[0]):
            logging.info('Calculating screen images...')
            residuals = inresiduals[:, :, sindx, newaxis].transpose([0, 2, 1]).reshape(N_piercepoints, N_times)
            mpm = multiprocManager(ncpu, _calculate_screen)
            for k in range(N_times):
                mpm.put([inscreen[:, k, sindx], residuals[:, k], pp,
                    N_piercepoints, k, east, north, up, T, Nx, Ny, sindx, height,
                    beta_val, r_0, is_phase])
            mpm.wait()
            for (k, ft, scr, xa, ya) in mpm.get():
                screen[:, :, k] = scr
                fitted_phase1[:, k] = ft
                if is_image_plane:
                    x[k, :] = xa
                    y[k, :] = ya
                else:
                    x[k, :] = xa - np.amin(xa) # remove offsets for each time slot
                    y[k, :] = ya - np.amin(ya)

            # Normalize piercepoint locations to extent calculated above
            if not is_image_plane:
                x *= extent[0]
                y *= extent[1]

            if min_val is None:
                vmin = np.min([np.amin(screen), np.amin(fitted_phase1)])
            else:
                vmin = min_val
            if max_val is None:
                vmax = np.max([np.amax(screen), np.amax(fitted_phase1)])
            else:
                vmax = max_val

            logging.info('Plotting screens...')
            mpm = multiprocManager(ncpu, _plot_frame)
            for k in range(N_times):
                mpm.put([screen[:, :, k], fitted_phase1[:, k], residuals[:, k],
                weights[:, k, sindx], x[k, :], y[k, :], k, lower, upper, vmin, vmax,
                source_names, show_source_names, station_names, sindx, root_dir,
                prestr, is_image_plane, midRA, midDec, station_order[0, k, sindx], is_phase])
            mpm.wait()
    else:
        logging.info('Calculating screen images...')
        residuals = inresiduals.transpose([0, 2, 1]).reshape(N_piercepoints, N_times)
        weights = weights.transpose([0, 2, 1]).reshape(N_piercepoints, N_times)
        mpm = multiprocManager(ncpu, _calculate_screen)
        for k in range(N_times):
            mpm.put([inscreen[:, k, :], residuals[:, k], pp[k, :, :],
                N_piercepoints, k, east, north, up, T, Nx, Ny, -1, height,
                beta_val, r_0, is_phase])
        mpm.wait()
        for (k, ft, scr, xa, ya) in mpm.get():
            screen[:, :, k] = scr
            fitted_phase1[:, k] = ft
            if is_image_plane:
                x[k, :] = xa
                y[k, :] = ya
            else:
                x[k, :] = xa - np.amin(xa) # remove offsets for each time slot
                y[k, :] = ya - np.amin(ya)

        # Normalize piercepoint locations to extent calculated above
        if not is_image_plane:
            x *= extent[0]
            y *= extent[1]

        if min_val is None:
            vmin = np.min([np.amin(screen), np.amin(fitted_phase1)])
        else:
            vmin = min_val
        if max_val is None:
            vmax = np.max([np.amax(screen), np.amax(fitted_phase1)])
        else:
            vmax = max_val

        logging.info('Plotting screens...')
        if show_source_names:
            pp_names = []
            for stat in station_names:
                for src in source_names:
                    pp_names.append('{0}_{1}'.format(src, stat))
        else:
            pp_names = source_names
        mpm = multiprocManager(ncpu, _plot_frame)
        for k in range(N_times):
            order = station_order[0]
            mpm.put([screen[:, :, k], fitted_phase1[:, k], residuals[:, k],
            weights[:, k], x[k, :], y[k, :], k, lower, upper, vmin, vmax,
            pp_names, show_source_names, station_names, -1, root_dir,
            prestr, is_image_plane, midRA, midDec, order, is_phase])
        mpm.wait()


def _fitPLaneLTSQ(XYZ):
    """
    Fits a plane to an XYZ point cloud

    Returns (a, b, c), where Z = aX + bY + c

    Parameters
    ----------
    XYZ : array
        point cloud

    Returns
    -------
    (a, b, c) : floats
        plane parameters, where Z = aX + bY + c

    """
    import numpy as np
    [rows, cols] = XYZ.shape
    G = np.ones((rows, 3))
    G[:, 0] = XYZ[:, 0]  #X
    G[:, 1] = XYZ[:, 1]  #Y
    Z = XYZ[:, 2]
    (a, b, c), resid, rank, s = np.linalg.lstsq(G, Z)
    return (a, b, c)


def run(soltab, resSoltab='', minZ=-3.2, maxZ=3.2, prefix='', remove_gradient=False,
    show_source_names=False, ncpu=0):
    """
    Plot screens (one plot is made per time and per station)

    Parameters
    ----------
    soltab : solution table
        Soltab containing screen.
    resSoltab : solution table, optional
        Soltab containing the screen residuals.
    minZ : float, optional
        Minimum value of colorbar scale.
    maxZ : float, optional
        Max value of colorbar scale.
    prefix : str, optional
        String to prepend to output plots.
    remove_gradient : bool, optional
        If True, remove gradient before plotting.
    show_source_names : bool, optional
        If True, source names are overplotted.
    ncpu : int, optional
        Number of CPUs to use. If 0, all are used.

    """
    import os
    import numpy as np

    screen_type = soltab.getType()
    if screen_type == 'phasescreen':
        is_phase = True
    else:
        is_phase = False
    logging.info('Using input {0} soltab {1}'.format(screen_type, soltab.name))

    # Get values from soltabs
    solset = soltab.getSolset()
    if resSoltab is '':
        try:
            # Look for residual soltab assuming standard naming conventions
            ressoltab = solset.getSoltab(soltab.name+'resid')
        except:
            logging.error('Could not find the soltab with associated screen residuals. '
                'Please specify it with the "resSoltab" argument.')
            return 1
    else:
        ressoltab = solset.getSoltab(resSoltab)
    logging.info('Using input screen residual soltab: {}'.format(ressoltab.name))
    screen = np.array(soltab.val)
    weights = np.array(soltab.weight)
    residuals = np.array(ressoltab.val)
    times = np.array(soltab.time)
    orders = np.array(ressoltab.weight)
    axis_names = soltab.getAxesNames()
    if len(screen.shape) > 3:
        # remove degenerate freq and pol axes
        if 'freq' in axis_names:
            freq_ind = axis_names.index('freq')
            screen = np.squeeze(screen, axis=freq_ind)
            weights = np.squeeze(weights, axis=freq_ind)
            residuals = np.squeeze(residuals, axis=freq_ind)
            orders = np.squeeze(orders, axis=freq_ind)
            axis_names.remove('freq')
        if 'pol' in axis_names:
            pol_ind = axis_names.index('pol')
            screen = np.squeeze(screen, axis=pol_ind)
            weights = np.squeeze(weights, axis=pol_ind)
            residuals = np.squeeze(residuals, axis=pol_ind)
            orders = np.squeeze(orders, axis=pol_ind)
            axis_names.remove('pol')

    # Rearrange to get order [dir, time, ant]
    dir_ind = axis_names.index('dir')
    time_ind = axis_names.index('time')
    ant_ind = axis_names.index('ant')
    screen = screen.transpose([dir_ind, time_ind, ant_ind])
    weights = weights.transpose([dir_ind, time_ind, ant_ind])
    residuals = residuals.transpose([dir_ind, time_ind, ant_ind])
    orders = orders.transpose([dir_ind, time_ind, ant_ind])

    # Collect station and source names and positions and times, making sure
    # that they are ordered correctly.
    source_names = soltab.dir[:]
    source_dict = solset.getSou()
    source_positions = []
    for source in source_names:
        source_positions.append(source_dict[source])
    station_names = soltab.ant
    station_dict = solset.getAnt()
    station_positions = []
    for station in station_names:
        station_positions.append(station_dict[station])
    height = soltab.obj._v_attrs['height']
    beta_val = soltab.obj._v_attrs['beta']
    r_0 = soltab.obj._v_attrs['r_0']
    pp = soltab.obj.piercepoint[:]
    if height == 0.0:
        midRA = soltab.obj._v_attrs['midra']
        midDec = soltab.obj._v_attrs['middec']
    else:
        midRA = 0.0
        midDec = 0.0

    if (minZ == 0 and maxZ == 0):
        min_val = None
        max_val = None
    else:
        min_val = minZ
        max_val = maxZ

    _make_screen_plots(pp, screen, residuals, weights, np.array(station_names),
        np.array(station_positions), np.array(source_names), times,
        height, orders, beta_val, r_0, prefix=prefix,
        remove_gradient=remove_gradient, show_source_names=show_source_names, min_val=min_val,
        max_val=max_val, is_phase=is_phase, midRA=midRA, midDec=midDec, ncpu=ncpu)

    return 0
