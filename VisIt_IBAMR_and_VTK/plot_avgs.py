#! /usr/bin/env python

'''This module is for loading and plotting the averaged magnitudes of each
velocity component in the flow past cylinders data.

Author: Christopher Strickland
Date: 6/16/2016
'''

# Compare different tower heights
# Compare different Re

# More complicated:
#   Compare the averaged x-component of the velocity vs. height to
#   the analytical solution, solving for the best fit value of the permeability
#       (the parameter alpha?). Use least squares.

from math import exp, log
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
plt.rcParams['image.cmap'] = 'viridis'
cmap = cm.get_cmap('viridis')

##################
# Get mu from input3d. It should be constant between runs
# dpdx - calculate derivative of pressure field (P) in VisIt at the last time
#   then average it over the entire domain
# length of cylinder is file_name_length/64
# height of domain = 1
# U is V in input3d (always 0.1)


def flow_model(z_ary,alpha,a,b,U,dpdx,mu):
    '''Evaluate the analytical 2D flow model at and above the porous tower.
    Flow is assumed to be moving in the x direction with constant velocity in
    that direction for a given value of z.

    Arguments:
        z_ary: position in the z direction. ndarray
            The base of the domain is -a, the top b, and the top of the cylinder
            is at 0.
        alpha: porous constant (to be fit)
        a: height of cylinder
        b: length of domain above cylinder
        U: u(b). This is v in input3d
        dpdx: dp/dx change in momentum constant (from data parsing)
        mu: constant, dynamic viscosity

    Returns:
        u(z): velocity of flow in the x direction at z
        '''

    # Calculate C and D constants and then get A and B based on these

    C = (dpdx*(-0.5*alpha**2*b**2+exp(log(alpha*b)-alpha*a)-exp(-alpha*a)+1) +
        U*alpha**2*mu)/(alpha**2*mu*(exp(log(alpha*b)-2*alpha*a)+alpha*b-
        exp(-2*alpha*a)+1))

    D = (dpdx*(exp(log(0.5*alpha**2*b**2)-2*alpha*a)+exp(log(alpha*b)-alpha*a)+
        exp(-alpha*a)-exp(-2*alpha*a)) - 
        exp(log(U*alpha**2*mu)-2*alpha*a))/(alpha**2*mu*
        (exp(log(alpha*b)-2*alpha*a)+alpha*b-exp(-2*alpha*a)+1))

    A = alpha*C - alpha*D
    B = C + D - dpdx/(alpha**2*mu)

    # Form solution array
    sol = np.zeros_like(z_ary)

    for n,z in enumerate(z_ary):
        if z > 0:
            #Region I
            sol[n] = z**2*dpdx/(2*mu) + A*z + B
        else:
            #Region 2
            if C > 0 and D > 0:
                sol[n] = exp(log(C)+alpha*z) + exp(log(D)-alpha*z) - dpdx/(alpha**2*mu)
            elif C <= 0 and D > 0:
                sol[n] = exp(log(D)-alpha*z) - dpdx/(alpha**2*mu)
            elif C > 0 and D <= 0:
                sol[n] = exp(log(C)+alpha*z) - dpdx/(alpha**2*mu)
            else:
                sol[n] = -dpdx/(alpha**2*mu)

    return sol


def resid(alpha,a,b,v,dpdx,mu,z_mesh,x_abs_avgs):
    '''This is mostly a wrapper around flow_model for the
    least squares algorithm'''

    # The analytical model expects the bottom of the domain to be at -a, where
    #   a is the height of the porous region. Thus we need to shift z_mesh,
    #   since it typically starts somewhere else (-0.5).
    z_mesh_model = z_mesh - z_mesh[0] - a
    # Get analytical u(z) values
    u = flow_model(z_mesh_model,alpha,a,b,v,dpdx,mu)

    # Return residuals
    return u - x_abs_avgs


def get_data(tower, re, length):
    '''Return averaged data from file using the established naming scheme'''
    if tower == 1:
        filename = 'MacrophyteAvgs/viz_IB3d_{}tower_Re{}_len{}_avgs.txt'.format(
            tower,re,length)
    else:
        filename = 'MacrophyteAvgs/viz_IB3d_{}towers_Re{}_len{}_avgs.txt'.format(
            tower,re,length)
    z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs = np.loadtxt(filename)
    with open(filename[:-8]+'dpdx.txt') as fobj:
        dpdx_avg = float(fobj.readline())
    return (z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs,dpdx_avg)



def get_data_spacing(spacing, re):
    '''Return averaged data from file using naming scheme for spacing runs.

    Arguments:
        spacing: 1, 2, 4, or 8 for 1/1, 1/2, 1/4, and 1/8 respectively
            (8 looks for original data, 1 tower, length 10/64)
        re: Reynolds number
    '''

    if spacing != 1 and spacing != 8:
        filename = 'MacrophyteAvgs/viz_IB3d_1tower_Re{}_1_{}_spacing_avgs.txt'.format(re,spacing)
    elif spacing == 1:
        filename = 'MacrophyteAvgs/viz_IB3d_1tower_Re{}_1_spacing_avgs.txt'.format(re)
    else:
        filename = 'macrophyteAvgs/viz_IB3d_1tower_Re{}_len10_avgs.txt'.format(re)
    z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs = np.loadtxt(filename)
    with open(filename[:-8]+'dpdx.txt') as fobj:
        dpdx_avg = float(fobj.readline())
    return (z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs,dpdx_avg)



def compare_heights(tower_num, re, length_list):
    '''Create a plot comparing averaged data across different tower heights.

    Arguments:
        tower_num: number of towers (int)
        re: Reynolds number (int)
        length_list: list of tower heights to compare across
    '''

    data_list = []
    for height in length_list:
        data_list.append(get_data(tower_num, re, height))

    # plot setup
    f, axarr = plt.subplots(3, sharex=True, figsize=(9,5.5))
    axarr[0].set_title(
        'Planar avgs. for fluid velocity around {} tower at Re {}.'.format(tower_num,re))
    axarr[0].set_ylabel('Avg. fluid speed in\n direction of flow, X(z)')
    axarr[0].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[0].set_ylim(0,0.01)
    axarr[1].set_ylabel(r'$\frac{\mathrm{Avg.} Y\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[1].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_ylabel(r'$\frac{\mathrm{Avg.} Z\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[2].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_xlabel('Z intercept of plane')

    # color setup
    color_list = np.linspace(0.85,0.05,len(length_list))

    for n,height in enumerate(length_list):
        for ii in range(3):
            if ii == 0:
                # X direction
                axarr[ii].plot(data_list[n][0],data_list[n][1],
                    label='Tower height: {}/64'.format(height),c=cmap(color_list[n]))
            else:
                # normalize Y and Z direction
                this_data = np.ma.array(data_list[n][ii+1]/data_list[n][1])
                masked_data = np.ma.masked_where(data_list[n][ii+1] < 1e-10, this_data)
                axarr[ii].plot(data_list[n][0],this_data,
                    label='Tower height: {}/64'.format(height),c=cmap(color_list[n]))
            # plot a vertical line at tower height
            top = height/64 - 0.5
            axarr[ii].axvline(x=top,color=cmap(color_list[n]),ls='--')
    for ii in range(3):
        # print legend
        leg = axarr[ii].legend(loc="upper right",fontsize=11)
        leg.get_frame().set_alpha(0.65)
    xlabels_float = axarr[2].get_xticks().tolist()
    xlabels_float -= data_list[0][0][0]
    axarr[2].set_xticklabels([str(item) for item in xlabels_float])
    plt.tight_layout()
    plt.show()



def compare_re(tower_num, re_list, length):
    '''Create a plot comparing averaged data across different Re numbers.

    Arguments:
        tower_num: number of towers (int)
        re: list of Reynolds numbers to compare across (ints)
        length: tower height (int)
    '''

    data_list = []
    for re in re_list:
        data_list.append(get_data(tower_num, re, length))

    # plot setup
    f, axarr = plt.subplots(3, sharex=True, figsize=(7.5,5.5))
    axarr[0].set_title(
        'Planar avgs. for {} tower fluid velocity by Re. Tower height: {}/64'.format(tower_num,length))
    axarr[0].set_ylabel('Avg. fluid speed in\n direction of flow, X(z)')
    axarr[0].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[1].set_ylabel(r'$\frac{\mathrm{Avg.} Y\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[1].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_ylabel(r'$\frac{\mathrm{Avg.} Z\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[2].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_xlabel('Z intercept of plane')

    # color setup
    color_list = np.linspace(0.85,0.05,len(re_list))
    cmap = cm.get_cmap('viridis')

    for n,re in enumerate(re_list):
        for ii in range(3):
            if ii == 0:
                # X direction
                axarr[ii].plot(data_list[n][0],data_list[n][ii+1],
                    label='Re = {}'.format(re),c=cmap(color_list[n]))
            else:
                # normalize Y and Z direction
                this_data = np.ma.array(data_list[n][ii+1]/data_list[n][1])
                masked_data = np.ma.masked_where(data_list[n][ii+1] < 1e-10, this_data)
                axarr[ii].plot(data_list[n][0],this_data,
                    label='Re = {}'.format(re),c=cmap(color_list[n]))
    for ii in range(3):
        # plot a vertical line at tower height
        axarr[ii].axvline(x=length/64 - 0.5,color='k',ls='--')
        # label it
        axarr[ii].text(length/64+0.035, 0.78, 'top of\ncylinder',
                       transform=axarr[ii].transAxes)
        # print legend
        if ii > 0:
            leg = axarr[ii].legend(loc="upper right",ncol=2,fontsize=11)
        else:
            leg = axarr[ii].legend(loc="lower right",ncol=2,fontsize=11)
        leg.get_frame().set_alpha(0.65)
    xlabels_float = axarr[2].get_xticks().tolist()
    xlabels_float -= data_list[0][0][0]
    axarr[2].set_xticklabels([str(item) for item in xlabels_float])
    plt.tight_layout()
    plt.show()



def compare_towernum(tower_list, re):
    '''Create a plot comparing different numbers of towers. Tower length = 10/64
    
    Arguments:
        tower_list: list of tower numbers (ints)
        re: Reynolds number'''

    data_list = []
    for tower_num in tower_list:
        data_list.append(get_data(tower_num, re, 10))

    # plot setup
    f, axarr = plt.subplots(3, sharex=True, figsize=(7.5,5.5))
    axarr[0].set_title(
        'Planar avgs. for fluid velocity by array size at Re {}.'.format(re))
    axarr[0].set_ylabel('Avg. fluid speed in\n direction of flow, X(z)')
    axarr[0].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[0].set_ylim(0,0.01)
    axarr[1].set_ylabel(r'$\frac{\mathrm{Avg.} Y\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[1].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_ylabel(r'$\frac{\mathrm{Avg.} Z\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[2].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_xlabel('Z intercept of plane')

    # color setup
    color_list = np.linspace(0.85,0.05,len(tower_list))
    cmap = cm.get_cmap('viridis')

    for n,num in enumerate(tower_list):
        for ii in range(3):
            if ii == 0:
                # X direction
                axarr[ii].plot(data_list[n][0],data_list[n][ii+1],
                    label='{} tower(s)'.format(num),c=cmap(color_list[n]))
            else:
                # normalize Y and Z direction
                this_data = np.ma.array(data_list[n][ii+1]/data_list[n][1])
                masked_data = np.ma.masked_where(data_list[n][ii+1] < 1e-10, this_data)
                axarr[ii].plot(data_list[n][0],this_data,
                    label='{} tower(s)'.format(num),c=cmap(color_list[n]))
    for ii in range(3):
        # plot a vertical line at tower height
        axarr[ii].axvline(x=10/64 - 0.5,color='k',ls='--')
        # label it
        axarr[ii].text(10/64+0.035, 0.78, 'top of\ncylinder',
                       transform=axarr[ii].transAxes)
        # print legend
        if ii > 0:
            leg = axarr[ii].legend(loc="upper right",ncol=2,fontsize=11)
        else:
            leg = axarr[ii].legend(loc="upper right",ncol=2,fontsize=11)
        leg.get_frame().set_alpha(0.65)
    xlabels_float = axarr[2].get_xticks().tolist()
    xlabels_float -= data_list[0][0][0]
    axarr[2].set_xticklabels([str(item) for item in xlabels_float])
    plt.tight_layout()
    plt.show()
    


def compare_all_Xdir(tower_list=[1,4,16], re_list=[0.2,1], length_list=[10]):
    '''Different plot for each tower length. Plots X-direction velocity only.'''

    plt.figure()
    N = len(re_list)*len(tower_list)
    axarr = []

    for ii, length in enumerate(length_list):
        axarr.append(plt.subplot(len(length_list),1,ii+1))
        axarr[ii].set_title('Avg. fluid velocity with tower length {}/64'.format(length))

        # Get color cycler
        cmap = plt.get_cmap('Set2')
        axarr[ii].set_color_cycle([cmap(k) for k in np.linspace(0,1,N)])
        color_cycle = axarr[ii]._get_lines.prop_cycler

        for jj, re in enumerate(re_list):
            for kk, tower in enumerate(tower_list):
                clrdict = next(color_cycle)
                # Get data
                z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs,dpdx_avg = \
                    get_data(tower,re,length)
                # plot data
                axarr[ii].plot(z_mesh,x_abs_avgs,color=clrdict['color'],
                    label='{} tower(s), Re {}'.format(tower,re))
        leg = axarr[ii].legend(loc="upper left",fontsize=12)
        leg.get_frame().set_alpha(0.65)

    # labels
    if len(length_list) > 3:
        axarr[2].set_ylabel('Average fluid velocity in direction of flow',fontsize=12)
    else:
        axarr[0].set_ylabel('Average fluid velocity in direction of flow',fontsize=12)
    axarr[-1].set_xlabel('Z intercept of plane')

    plt.show()



def compare_spacing(spacing_list, re):
    '''Create a plot comparing tower spacing and Reynolds number

    Arguments:
        spacing_list: list of spacings (1,2,4,8) for (1/1,1/2,1/4,1/8)
        re: Reynolds number
    '''

    data_list = []
    for spacing in spacing_list:
        data_list.append(get_data_spacing(spacing, re))

    # plot setup
    f, axarr = plt.subplots(3, sharex=True, figsize=(9,5.5))
    axarr[0].set_title(
        'Planar avgs. for tower fluid velocity by spacing. Re = {}, Tower height = 10/64'.format(re))
    axarr[0].set_ylabel('Avg. fluid speed in\n direction of flow, X(z)')
    axarr[0].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[1].set_ylabel(r'$\frac{\mathrm{Avg.} Y\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[1].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_ylabel(r'$\frac{\mathrm{Avg.} Z\ \mathrm{flow\ speed}}{\mathrm{Avg.} X\ \mathrm{flow\ speed}}$',
        fontsize=14)
    axarr[2].set_xlim(data_list[0][0][0],data_list[0][0][-1])
    axarr[2].set_xlabel('Z intercept of plane')

    # color setup
    color_list = np.linspace(0.85,0.05,len(spacing_list))
    cmap = cm.get_cmap('viridis')

    for n,spacing in enumerate(spacing_list):
        for ii in range(3):
            if ii == 0:
                # X direction
                axarr[ii].plot(data_list[n][0],data_list[n][ii+1],
                    label='Spacing = 1/{}'.format(spacing),c=cmap(color_list[n]))
            else:
                # normalize Y and Z direction
                this_data = np.ma.array(data_list[n][ii+1]/data_list[n][1])
                masked_data = np.ma.masked_where(data_list[n][ii+1] < 1e-10, this_data)
                axarr[ii].plot(data_list[n][0],this_data,
                    label='Spacing = 1/{}'.format(spacing),c=cmap(color_list[n]))
    for ii in range(3):
        # plot a vertical line at tower height
        axarr[ii].axvline(x=10/64 - 0.5,color='k',ls='--')
        # label it
        axarr[ii].text(10/64+0.035, 0.78, 'top of\ncylinder',
                       transform=axarr[ii].transAxes)
        # print legend
        if ii > 0:
            leg = axarr[ii].legend(loc="upper right",ncol=2,fontsize=11)
        else:
            leg = axarr[ii].legend(loc="lower right",ncol=2,fontsize=11)
        leg.get_frame().set_alpha(0.65)
    xlabels_float = axarr[2].get_xticks().tolist()
    xlabels_float -= data_list[0][0][0]
    axarr[2].set_xticklabels([str(item) for item in xlabels_float])
    plt.tight_layout()
    plt.show()



def fit_model(tower_num, re, length):
    '''Find the alpha (permeability) that provides the best fit of the analytical
    model to the data.

    Arguments:
        tower_num: number of towers
        re: Reynolds number
        length: tower height'''

    # Parse the Reynolds number into component parameters
    rho = 1000 # density of water is 1000 kg/m**3
    v = 0.1 # max velocity of fluid
    L = 1 # characteristic length
    mu = rho*v*L/re

    z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs,dpdx_avg = get_data(tower_num,re,length)

    a = length/64
    b = L - a

    # Solve nonlinear least squares
    result_obj = least_squares(resid,1.,bounds=(0,np.inf),
                               args=(a,b,v,dpdx_avg,mu,z_mesh,x_abs_avgs))

    # Get fitted analytical solution
    # The analytical model expects the bottom of the domain to be at -a, where
    #   a is the height of the porous region. Thus we need to shift z_mesh,
    #   since it typically starts somewhere else (-0.5).
    z_mesh_model = z_mesh - z_mesh[0] - a

    model = flow_model(z_mesh_model,float(result_obj.x),a,b,v,dpdx_avg,mu)

    #import pdb; pdb.set_trace()

    # Plot optimal analytical solution with data
    plt.figure()
    ax = plt.axes()
    plt.title('Flow velocity: data vs. optimal model')
    plt.hold(True)
    plt.plot(z_mesh,model,'r--',label='best-fit Brinkmann model')
    plt.scatter(z_mesh,x_abs_avgs,label='x-velocity data')
    plt.hold(False)
    plt.text(0.05,0.88,r'Optimal $\alpha$ = {}'.format(float(result_obj.x))+
        '\nNumber of towers = {}\nRe = {}\n'.format(tower_num,re)+
        'Tower length = {}'.format(length),
        ha='left',va='center',transform=ax.transAxes,fontsize=14)
    plt.xlabel('Z')
    plt.ylabel('velocity')
    plt.legend(loc='lower right')
    plt.show()



def fit_model_loop(tower_num,re, length):
    '''Find the alpha (permeability) that provides the best fit of the analytical
    model to the data. At least one parameter must be a list to be iterated over.

    Arguments:
        tower_num: number of towers
        re: Reynolds number
        length: tower height'''

    # Make all arguments lists
    if not isinstance(tower_num,list):
       tower_num = [tower_num]
    if not isinstance(re,list):
        re = [re]
    if not isinstance(length,list):
        length = [length]

    # Parse the Reynolds number into component parameters
    rho = 1000 # density of water is 1000 kg/m**3
    v = 0.1 # max velocity of fluid
    L = 1 # characteristic length

    result_obj = []
    model = []

    plt.figure()
    ax = plt.axes()
    plt.hold(True)
    plt.title('Flow velocity: simulation data vs. optimal model',fontsize="x-large")

    # Get color cycler
    N = len(tower_num)*len(re)*len(length)
    cmap = plt.get_cmap('Set2')
    ax.set_color_cycle([cmap(k) for k in np.linspace(0,1,N)])
    color_cycle = ax._get_lines.prop_cycler

    # zoomed view
    axins = inset_axes(ax,"50%","30%",loc=6)

    # Iterate over the lists
    for ii,tower in enumerate(tower_num):
        for jj,r in enumerate(re):
            for kk,l in enumerate(length):

                mu = rho*v*L/r

                # Get data
                z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs,dpdx_avg = \
                    get_data(tower,r,l)

                a = l/64
                b = L - a

                # Solve nonlinear least squares, just fitting up to tower height
                zind = np.argwhere(z_mesh>a-0.5)[0].flatten()
                result_obj = least_squares(resid,1.,bounds=(0,np.inf),
                                    args=(a,b,v,dpdx_avg,mu,
                                    z_mesh[:zind],x_abs_avgs[:zind]))

                # Get fitted analytical solution
                # The analytical model expects the bottom of the domain to be at -a, where
                #   a is the height of the porous region. Thus we need to shift z_mesh,
                #   since it typically starts somewhere else (-0.5).
                z_mesh_model = z_mesh - z_mesh[0] - a

                model = flow_model(
                    z_mesh_model,float(result_obj.x),a,b,v,dpdx_avg,mu)

                # Plot optimal analytical solution with data
                line, = ax.plot(z_mesh,model,'k--',label='best-fit Brinkmann model')
                axins.plot(z_mesh,model,'k--',label='best-fit Brinkmann model')
                if ii+jj+kk > 0:
                    line.set_label('_no label')
                clrdict = next(color_cycle)
                ax.scatter(z_mesh,x_abs_avgs,
                    label='{} tower, Re {}, len {}/64. '.format(tower,r,l)+\
                    r'$\alpha$ = '+'{:.2f}'.format(result_obj.x[0]),
                    color=clrdict['color'],alpha=0.65)
                axins.scatter(z_mesh,x_abs_avgs,
                    label='{} tower, Re {}, len {}/64. '.format(tower,r,l)+\
                    r'$\alpha$ = '+'{:.2f}'.format(result_obj.x[0]),
                    color=clrdict['color'],alpha=0.65)
                # Plot where the top of the cylinder is
                zind = np.argwhere(z_mesh<=a-0.5)
                ax.plot(z_mesh[zind[-1]],x_abs_avgs[zind[-1]],'*',
                    color=clrdict['color'],mec=clrdict['color'],ms=16,mew=2)
                axins.plot(z_mesh[zind[-1]],x_abs_avgs[zind[-1]],'*',
                    color=clrdict['color'],mec=clrdict['color'],ms=16,mew=2)
                # Print out optimal alpha values
                print(r'{} towers, Re {}, len {}, $\alpha$ = {}'.format(
                    tower,r,l,result_obj.x
                ))

    ax.set_xlabel('Z intercept of plane',fontsize='large')
    ax.set_ylabel('Velocity in direction of flow',fontsize='large')
    ax.set_xlim([-0.5,0.425])
    ax.set_ylim(ymin=0)
    #plt.ylim([0,0.008])
    leg = ax.legend(loc='upper left',fontsize='medium')
    leg.get_frame().set_alpha(0.65)

    axins.set_xlim(-0.5,0.15)
    axins.set_ylim(0,0.005)
    plt.yticks(visible=False)
    mark_inset(ax,axins,loc1=1,loc2=3,fc="none",ec="0.45")
    plt.show()



def journal_figures():
    '''Plot and save journal figures'''

    '''Fit model for variable tower length'''
    tower_num = 1 # one tower
    re = 1 # Re number
    length = [10,20,30,40] # lengths

    # Make all arguments lists
    if not isinstance(tower_num,list):
       tower_num = [tower_num]
    if not isinstance(re,list):
        re = [re]
    if not isinstance(length,list):
        length = [length]

    # Parse the Reynolds number into component parameters
    rho = 1000 # density of water is 1000 kg/m**3
    v = 0.1 # max velocity of fluid
    L = 1 # characteristic length

    result_obj = []
    model = []

    plt.figure()
    ax = plt.axes()
    plt.hold(True)
    plt.title('Flow velocity: simulation data vs. optimal model',fontsize="x-large")

    # Get color cycler
    N = len(tower_num)*len(re)*len(length)
    cmap = plt.get_cmap('Set2')
    ax.set_color_cycle([cmap(k) for k in np.linspace(0,1,N)])
    color_cycle = ax._get_lines.prop_cycler

    # zoomed view
    axins = inset_axes(ax,"50%","30%",loc=6)

    # Iterate over the lists
    for ii,tower in enumerate(tower_num):
        for jj,r in enumerate(re):
            for kk,l in enumerate(length):

                mu = rho*v*L/r

                # Get data
                z_mesh,x_abs_avgs,y_abs_avgs,z_abs_avgs,dpdx_avg = \
                    get_data(tower,r,l)

                a = l/64
                b = L - a

                # Solve nonlinear least squares, just fitting up to tower height
                zind = np.argwhere(z_mesh>a-0.5)[0].flatten()
                result_obj = least_squares(resid,1.,bounds=(0,np.inf),
                                    args=(a,b,v,dpdx_avg,mu,
                                    z_mesh[:zind],x_abs_avgs[:zind]))

                # Get fitted analytical solution
                # The analytical model expects the bottom of the domain to be at -a, where
                #   a is the height of the porous region. Thus we need to shift z_mesh,
                #   since it typically starts somewhere else (-0.5).
                z_mesh_model = z_mesh - z_mesh[0] - a

                model = flow_model(
                    z_mesh_model,float(result_obj.x),a,b,v,dpdx_avg,mu)

                # Plot optimal analytical solution with data
                line, = ax.plot(z_mesh,model,'k--',label='best-fit Brinkmann model')
                axins.plot(z_mesh,model,'k--',label='best-fit Brinkmann model')
                if ii+jj+kk > 0:
                    line.set_label('_no label')
                clrdict = next(color_cycle)
                ax.scatter(z_mesh,x_abs_avgs,
                    label='{} tower, Re {}, len {}/64. '.format(tower,r,l)+\
                    r'$\alpha$ = '+'{:.2f}'.format(result_obj.x[0]),
                    color=clrdict['color'],alpha=0.65)
                axins.scatter(z_mesh,x_abs_avgs,
                    label='{} tower, Re {}, len {}/64. '.format(tower,r,l)+\
                    r'$\alpha$ = '+'{:.2f}'.format(result_obj.x[0]),
                    color=clrdict['color'],alpha=0.65)
                # Plot where the top of the cylinder is
                zind = np.argwhere(z_mesh<=a-0.5)
                ax.plot(z_mesh[zind[-1]],x_abs_avgs[zind[-1]],'*',
                    color=clrdict['color'],mec=clrdict['color'],ms=16,mew=2)
                axins.plot(z_mesh[zind[-1]],x_abs_avgs[zind[-1]],'*',
                    color=clrdict['color'],mec=clrdict['color'],ms=16,mew=2)
                # Print out optimal alpha values
                print(r'{} towers, Re {}, len {}, $\alpha$ = {}'.format(
                    tower,r,l,result_obj.x
                ))

    ax.set_xlabel('Z intercept of plane',fontsize='large')
    ax.set_ylabel('Velocity in direction of flow',fontsize='large')
    ax.set_xlim([-0.5,0.425])
    ax.set_ylim(ymin=0)
    #plt.ylim([0,0.008])
    leg = ax.legend(loc='upper left',fontsize='medium')
    leg.get_frame().set_alpha(0.65)

    axins.set_xlim(-0.5,0.15)
    axins.set_ylim(0,0.005)
    plt.yticks(visible=False)
    mark_inset(ax,axins,loc1=1,loc2=3,fc="none",ec="0.45")
    plt.show()



def main():
    '''Interactive menu here'''
    pass

if __name__ == "__main__":
    main()