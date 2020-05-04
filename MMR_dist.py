import numpy as np
from astropy import constants as const
from astropy import units as u
import pandas as pd
import matplotlib.cm as cm
from pylab import *
from matplotlib import gridspec


csv_file = 'planets3.csv'

class MMR_dist():
    def __init__(self,csv_file=csv_file):
        self.df_data = pd.read_csv(csv_file)

        self.names = self.df_data['NAME']
        self.masses = self.df_data['MASS (Mjup)']
        self.periods = self.df_data['PER (days)']
        self.Rstars = self.df_data['RSTAR (Rsun)']
        self.Mstars = self.df_data['MSTAR (Msun)']
        self.a_semi_au = self.df_data['A (au)']
        self.a_semi_Rstar = self.df_data['A (au)'] * u.au.to(u.Rsun) / self.Rstars
        self.stars = self.df_data['STAR']
        self.ecc = self.df_data['ECC']

    def clean_dat(self, params):

        names_temp = np.copy(self.names)
        masses_temp = np.copy(self.masses)
        periods_temp = np.copy(self.periods)
        Rstars_temp = np.copy(self.Rstars)
        Mstars_temp = np.copy(self.Mstars)
        a_semi_au_temp = np.copy(self.a_semi_au)
        a_semi_Rstar_temp = np.copy(self.a_semi_Rstar)
        stars_temp = np.copy(self.stars)
        ecc_temp = np.copy(self.ecc)

        for n in range(len(params)):

            if params[n] in params:
                if params[n] == 'Planet Mass':
                    where_none = np.where(masses_temp != 10000)[0]

                if params[n] == 'Stellar Radius':
                    where_none = np.where(Rstars_temp != 10000)[0]

                if params[n] == 'Stellar Mass':
                    where_none = np.where(Mstars_temp != 10000)[0]

                if params[n] == 'Planet Eccentricity':
                    where_none = np.where(ecc_temp != 10000)[0]

                    # self.names = self.names[where_none1]
                    # self.masses = self.masses[where_none1]
                    # self.periods = self.periods[where_none1]
                    # self.Rstars = self.Rstars[where_none1]
                    # self.Mstars = self.Mstars[where_none1]
                    # self.a_semi = self.a_semi[where_none1]
                    # self.stars = self.stars[where_none1]
                    # self.ecc = self.ecc[where_none1]

                names_temp = names_temp[where_none]
                masses_temp = masses_temp[where_none]
                periods_temp = periods_temp[where_none]
                Rstars_temp = Rstars_temp[where_none]
                Mstars_temp = Mstars_temp[where_none]
                a_semi_au_temp = a_semi_au_temp[where_none]
                a_semi_Rstar_temp = a_semi_Rstar_temp[where_none]
                stars_temp = stars_temp[where_none]
                ecc_temp = ecc_temp[where_none]

        self.names = names_temp
        self.masses = masses_temp
        self.periods = periods_temp
        self.Rstars = Rstars_temp
        self.Mstars = Mstars_temp
        self.a_semi_au = a_semi_au_temp
        self.a_semi_Rstar = a_semi_Rstar_temp
        self.stars = stars_temp
        self.ecc = ecc_temp

        #return self.names, self.masses, self.periods, self.Rstars ,self.Mstars, self.a_semi, self.stars, self.ecc

    def find_multiples(self):
        names_mult = []
        masses_mult = []
        periods_mult = []
        a_semi_au_mult = []
        a_semi_Rstar_mult = []
        stars_mult = []
        ecc_mult = []
        Mstars_mult = []
        Rstars_mult = []

        used = ['none']
        for i in range(1, len(self.stars)):
            if self.stars[i] not in used:
                where_star = np.where(self.stars == self.stars[i])[0]

                if len(where_star) > 1:
                    names_mult.append(self.names[where_star])
                    masses_mult.append(self.masses[where_star])
                    periods_mult.append(self.periods[where_star])
                    a_semi_au_mult.append(self.a_semi_au[where_star])
                    a_semi_Rstar_mult.append(self.a_semi_Rstar[where_star])
                    ecc_mult.append(self.ecc[where_star])
                    Mstars_mult.append(self.Mstars[where_star])
                    Rstars_mult.append(self.Rstars[where_star])

                    stars_mult.append(self.stars[i])
                    used.append(self.stars[i])
                    continue

        names_mult_Pless = []
        masses_mult_Pless = []
        periods_mult_Pless = []
        a_semi_au_mult_Pless = []
        a_semi_Rstar_mult_Pless = []
        ecc_mult_Pless = []
        Mstars_mult_Pless = []
        Rstars_mult_Pless = []
        for j in range(len(periods_mult)):
            where_P_less = np.where(periods_mult[j] <= 150.)[0]

            if len(where_P_less) > 1:
                periods_insert, masses_insert, names_insert, a_semi_au_insert, a_semi_Rstar_insert, ecc_insert, Mstars_insert, Rstars_insert = zip(*sorted(zip(periods_mult[j], masses_mult[j], names_mult[j], a_semi_au_mult[j], a_semi_Rstar_mult[j], ecc_mult[j], Mstars_mult[j], Rstars_mult[j])))

                names_mult_Pless.append(names_insert)
                masses_mult_Pless.append(masses_insert)
                periods_mult_Pless.append(periods_insert)
                a_semi_au_mult_Pless.append(a_semi_au_insert)
                a_semi_Rstar_mult_Pless.append(a_semi_Rstar_insert)
                ecc_mult_Pless.append(ecc_insert)
                Mstars_mult_Pless.append(Mstars_insert)
                Rstars_mult_Pless.append(Rstars_insert)

        self.names_mult_Pless = names_mult_Pless
        self.masses_mult_Pless = masses_mult_Pless
        self.periods_mult_Pless = periods_mult_Pless
        self.a_semi_au_mult_Pless = a_semi_au_mult_Pless
        self.a_semi_Rstar_mult_Pless = a_semi_Rstar_mult_Pless
        self.ecc_mult_Pless = ecc_mult_Pless
        self.Mstars_mult_Pless = Mstars_mult_Pless
        self.Rstars_mult_Pless = Rstars_mult_Pless

    def calc_axes(self):
        period_ratio = []
        inner_periods = []
        inner_masses = []
        inner_a_semi_au = []
        outer_a_semi_au = []
        inner_a_semi_Rstar = []
        outer_a_semi_Rstar = []
        a_semi_ratio = []
        outer_masses = []
        mass_ratio = []
        ecc_ratio = []
        inner_ecc = []
        outer_ecc = []
        plot_mstars = []
        plot_rstars = []
        for k in range(len(self.periods_mult_Pless)):

            if len(self.periods_mult_Pless[k]) > 1:

                if np.float(self.masses_mult_Pless[k][0]) > 0:
                    p_ratio = np.float(self.periods_mult_Pless[k][1]) / np.float(self.periods_mult_Pless[k][0])
                    m_ratio = np.float(self.masses_mult_Pless[k][1]) / np.float(self.masses_mult_Pless[k][0])
                    a_ratio = np.float(self.a_semi_au_mult_Pless[k][1]) / np.float(self.a_semi_au_mult_Pless[k][0])

                    # if self.ecc_mult_Pless[k][1] == self.ecc_mult_Pless[k][0]:
                    #     e_ratio == 1
                    # else:
                    #     if self.ecc_mult_Pless[k][0] > 0:
                    #         e_ratio = np.float(self.ecc_mult_Pless[k][1]) / np.float(self.ecc_mult_Pless[k][0])
                    #     if self.ecc_mult_Pless[k][0] > 0:
                    #     e_ratio = np.float(self.ecc_mult_Pless[k][1]) / np.float(self.ecc_mult_Pless[k][0])

                    period_ratio.append(p_ratio)
                    inner_periods.append(self.periods_mult_Pless[k][0])
                    inner_masses.append(self.masses_mult_Pless[k][0])
                    outer_masses.append(self.masses_mult_Pless[k][1])
                    inner_a_semi_au.append(self.a_semi_au_mult_Pless[k][0])
                    outer_a_semi_au.append(self.a_semi_au_mult_Pless[k][1])
                    inner_a_semi_Rstar.append(self.a_semi_Rstar_mult_Pless[k][0])
                    outer_a_semi_Rstar.append(self.a_semi_Rstar_mult_Pless[k][1])
                    a_semi_ratio.append(a_ratio)
                    mass_ratio.append(m_ratio)
                    inner_ecc.append(self.ecc_mult_Pless[k][0])
                    outer_ecc.append(self.ecc_mult_Pless[k][1])
                    #ecc_ratio.append(e_ratio)
                    plot_mstars.append(self.Mstars_mult_Pless[k][0])
                    plot_rstars.append(self.Rstars_mult_Pless[k][0])

        self.period_ratio = period_ratio
        self.inner_periods = inner_periods
        self.inner_masses = inner_masses
        self.inner_a_semi_au = inner_a_semi_au
        self.outer_a_semi_au = outer_a_semi_au
        self.inner_a_semi_Rstar = inner_a_semi_Rstar
        self.outer_a_semi_Rstar = outer_a_semi_Rstar
        self.a_semi_ratio = a_semi_ratio
        self.outer_masses = outer_masses
        self.mass_ratio = mass_ratio
        self.ecc_ratio = ecc_ratio
        self.inner_ecc = inner_ecc
        self.outer_ecc = outer_ecc
        self.plot_mstars = plot_mstars
        self.plot_rstars = plot_rstars

    def plot_scatter(self, ymax=30, x_axis='period', x_unit='au', z_axis='m_ratio'):

        if x_axis == 'period':
            x = self.inner_periods
        if x_axis == 'a_semi':
            if x_unit == 'au':
                x = self.inner_a_semi_au
            if x_unit == 'Rstar':
                x = self.inner_a_semi_Rstar

        if z_axis == 'm_ratio':
            z = self.mass_ratio
            cbar_min = 0
            cbar_max = 30
        if z_axis == 'inner_ecc':
            z = self.inner_ecc
            cbar_min = 0
            cbar_max = max(z)
        if z_axis == 'outer_ecc':
            z = self.outer_ecc
            cbar_min = 0
            cbar_max = max(z)
        if z_axis == 'ecc_ratio':
            z = self.ecc_ratio
            cbar_min = 0
            cbar_max = 10
        if z_axis == 'a_ratio':
            z = self.a_semi_ratio
            cbar_min = 0
            cbar_max = max(z)
        if z_axis == 'stellar_mass':
            z = self.plot_mstars
            cbar_min = 0
            cbar_max = max(z)


        y = self.period_ratio


        fig = figure(num=None, figsize=(5, 5), facecolor='w', dpi=300)

        gs1 = gridspec.GridSpec(1, 1)
        gs1.update(left=0.0, right=1., hspace=0.0)
        ax = subplot(gs1[0, 0])

        font_size = 'medium'

        # scatter(x, y, c=z, cmap=cm.rainbow, s=np.pi * 0.5 ** 2)

        ax.set_title('cbar = ' + z_axis, fontsize=font_size, style='normal', family='sans-serif')
        ax.set_ylabel(r'P$_{2}$ / P$_{1}$', fontsize=font_size, style='normal', family='sans-serif')
        ax.set_ylim([0, ymax])
        if x_axis == 'period':
            ax.set_xlabel(r'P$_{1}$', fontsize=font_size, style='normal', family='sans-serif')
            ax.set_xlim([0, 25])
            ax.plot([0, 50], [1.5, 1.5], '--', color='#000000', lw=1, alpha=0.3)
        if x_axis == 'a_semi':
            if x_unit == 'au':
                ax.set_xlabel(r'a$_{1}$ (au)', fontsize=font_size, style='normal', family='sans-serif')
                ax.set_xlim([0, max(x)+0.025])
                ax.plot([0.00, max(x)+0.025], [1.5, 1.5], '--', color='#000000', lw=1, alpha=0.3)
            if x_unit == 'Rstar':
                ax.set_xlabel(r'a$_{1}$ (R$_{star}$)', fontsize=font_size, style='normal', family='sans-serif')
                ax.set_xlim([0, max(x)+3])
                ax.plot([0.00, max(x)+3], [1.5, 1.5], '--', color='#000000', lw=1, alpha=0.3)

        sc = ax.scatter(x, y, c=z, cmap=cm.rainbow, lw=0.85, s=np.pi * 1.2 ** 2, vmin=cbar_min, vmax=cbar_max)
        colorbar(sc)
        # cbar.ax.set_label(r'Mass of Inner Planet (M$_{J}$)', fontsize=font_size, style='normal', family='sans-serif')

        ax.tick_params(axis='both', labelsize=font_size, direction='in', top=True, right=True)
        if x_axis == 'period':
            fig.savefig('/Users/lbiddle/PycharmProjects/K2/Resonances/ratios_'+z_axis+'_period.pdf', format='pdf', bbox_inches='tight')
        if x_axis == 'a_semi':
            if x_unit == 'au':
                fig.savefig('/Users/lbiddle/PycharmProjects/K2/Resonances/ratios_'+z_axis+'_au.pdf', format='pdf', bbox_inches='tight')
            if x_unit == 'Rstar':
                fig.savefig('/Users/lbiddle/PycharmProjects/K2/Resonances/ratios_'+z_axis+'_Rstar.pdf', format='pdf', bbox_inches='tight')
        close()


MMR = MMR_dist()

clean_params = ['Planet Mass','Stellar Radius','Stellar Mass','Planet Eccentricity']
MMR.clean_dat(clean_params)

MMR.find_multiples()

MMR.calc_axes()

z_axes = ['m_ratio','inner_ecc','outer_ecc','stellar_mass','a_ratio']
x_axes = ['period','a_semi']
x_units = ['au','Rstar']

for g in z_axes:
    for h in x_axes:
        if h == 'a_semi':
            for p in x_units:
                MMR.plot_scatter(ymax=30, x_axis=h, x_unit=p, z_axis=g)
        if h == 'period':
            MMR.plot_scatter(ymax=30, x_axis=h, z_axis=g)




