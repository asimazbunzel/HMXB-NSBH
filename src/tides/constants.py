# Values used in this module were taken from MESA const module,
# located in $MESA_DIR/const/public/const.def
# Units are in CGS


# math constants
pi = 3.1415926535897932384626433832795028841971693993751e0
pi4 = 4*pi
eulercon = 0.577215664901532861e0
ln2 = 6.9314718055994529e-01 # = log_cr(2d0)
ln3 = 1.0986122886681096e+00 # = log_cr(3d0)
ln10 = 2.3025850929940455e+00 # = log_cr(10d0)
a2rad = pi/180.0e0 # angle to radians
rad2a = 180.0e0/pi # radians to angle
one_third = 1e0/3e0
two_thirds = 2e0/3e0
ln4pi3 = 1.4324119583011810e0 # = log_cr(4*pi/3)
two_13 = 1.2599210498948730e0 # = pow_cr(2d0,1d0/3d0)
four_13 = 1.5874010519681994e0 # = pow_cr(4d0,1d0/3d0)

standard_cgrav = 6.67428e-8 # gravitational constant (g^-1 cm^3 s^-2)
planck_h = 6.62606896e-27 # Planck's constant (erg s)
hbar = planck_h / (2*pi)
qe = 4.80320440e-10 # electron charge (esu == (g cm^3 s^-2)^(1/2))
avo = 6.02214179e23 # Avogadro's constant (mole^-1)
clight = 2.99792458e10 # speed of light in vacuum (cm s^-1)
kerg = 1.3806504e-16 # Boltzmann's constant (erg K^-1)
boltzm = kerg
cgas = boltzm*avo # ideal gas constant; erg/K
kev = 8.617385e-5 # converts temp to ev (ev K^-1)
amu = 1e0/avo# atomic mass unit (g)

mn = 1.6749286e-24 # neutron mass (g)
mp = 1.6726231e-24 # proton mass (g)
me = 9.1093826e-28 # (was 9.1093897d-28) electron mass (g)
rbohr = hbar*hbar/(me * qe * qe) # Bohr radius (cm)
fine = qe*qe/(hbar*clight) # fine structure constant
hion = 13.605698140e0 # hydrogen ionization energy (eV)
ev2erg = 1.602176487e-12 # electron volt (erg)
mev_to_ergs = 1e6*ev2erg
mev_amu = mev_to_ergs/amu
Qconv = mev_to_ergs*avo

boltz_sigma = 5.670400e-5 # boltzmann's sigma = a*c/4 (erg cm^-2 K^-4 s^-1)
crad = boltz_sigma*4/clight # = radiation density constant, a (erg cm^-3 K^-4);
                            # Prad = crad * T^4 / 3 # approx = 7.5657e-15

ssol = boltz_sigma
asol = crad
weinlam = planck_h*clight/(kerg * 4.965114232e0)
weinfre = 2.821439372e0*kerg/planck_h
rhonuc = 2.342e14 # density of nucleus (g cm^3)

# solar age, L, and R values from Bahcall et al, ApJ 618 (2005) 1049-1056.
msol = 1.9892e33  # solar mass (g)  <<< gravitational mass, not baryonic
rsol = 6.9598e10 # solar radius (cm)
lsol = 3.8418e33  # solar luminosity (erg s^-1)
agesol = 4.57e9  # solar age (years)
Msun = msol
Rsun = rsol
Lsun = lsol
Msun33 = msol*1e-33
Rsun11 = rsol*1e-11
Lsun33 = lsol*1e-33
ly = 9.460528e17 # light year (cm)
pc = 3.261633e0 * ly # parsec (cm)
secyer = 3.1558149984e7 # seconds per year
dayyer = 365.25e0 # days per year

Teffsol=5777.0e0
loggsol=4.4378893534131256e0
Teffsun=Teffsol
loggsun=loggsol

mbolsun = 4.74e0  # Bolometric magnitude of the Sun IAU Resolution B1 2015
mbolsol = mbolsun

m_earth = 5.9764e27 # earth mass (g) = 3.004424e-6 Msun
r_earth = 6.37e8 # earth radius (cm)
au = 1.495978921e13 # astronomical unit (cm)

m_jupiter = 1.8986e30 # jupiter mass (g) = 0.954454d-3 Msun
r_jupiter = 6.9911e9 # jupiter mean radius (cm)
semimajor_axis_jupiter = 7.7857e13 # jupiter semimajor axis (cm)
