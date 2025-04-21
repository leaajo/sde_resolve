# import des modules
import numpy as np  
import scipy.stats as stats  # Distributions statistiques
import scipy.integrate as spi  # Pour les équations différentielles ordinaires
import matplotlib.pyplot as plt  
from math import sqrt 

from scipy.stats import norm

def graphique(tableau_poids,temps_final, finesse_pas):
   
   nombre_observations = int(temps_final / finesse_pas) + 1
   temps = np.linspace(0, temps_final, nombre_observations)

    # graphique
   plt.figure(figsize=(12, 6))
   for j in range(N):
     plt.plot(temps, tableau_poids[:, j], label=f'Poisson {j+1}')

   plt.xlabel('Temps')
   plt.ylabel('Poids')
   plt.title('Évolution du poids des poissons au fil du temps')
   plt.legend()
   plt.show()

# FONCTIONS MODELE LOGI 1
def b_logistique1(x, k, N, x_max=25):
    return (0.05+0.1*(k/N))*x*(1-x/x_max)

def sigma_logistique1(c, x):
    return c * x

def sigma_der_logi1(c,x):
    return c

# FONCTIONS MODELE MULLER FEUGA
def b_muller(x, k, N, a = 0.001, b = 0.001):
    return (a*x+b)**(0.1+1.5*(k/N))

def sigma_muller(c, x, alpha = 1):
    return c * (x**alpha)

def sigma_der_muller(c, x, alpha = 1) :
    return alpha*c*x**(alpha-1)

# FONCTION MODELE LOGI 2 A
def b_logi2A(x, k, N, x_max=25):
    return (0.01+0.1*(k/N))*x*x*(1-x/x_max)

def sigma_logi2A(c, x):
    return c * x

def sigma_der_logi2A(c,x):
    return c

# FONCTION MODELE LOGI 2 B
def b_logi2B(x, k, N, x_max=25):
    return (0.01+0.1*(k/N))*x*x*((1-x/x_max)**2)

def sigma_logi2B(c, x, x_max=25):
    return c * x*(1-x/x_max)

def sigma_der_logi2B(c, x, x_max=25):
    return c*x*(1-2*x/x_max)

# METHODE EULER EXPLICITE
def euler_explicite(N, temps_final, finesse_pas, c, mu, modele):
    nombre_observations =int(temps_final / finesse_pas) + 1
    # Initialisation du tableau
    tableau_poids = np.zeros((nombre_observations, N))
    X_0 = np.random.uniform(0.5, 2, N)  
    tableau_poids[0, :] = X_0

    # Remplissage du tableau
    for i in range(1, nombre_observations):
    
        # On ordonne les poids
     poids_ordonnes = tableau_poids[i-1, :]
     poids_ordonnes = np.sort(poids_ordonnes)

     for j in range(N):
          #k = np.where(poids_ordonnes == tableau_poids[i-1, j])[0][0]
         k = np.searchsorted(poids_ordonnes, tableau_poids[i-1, j])
         x = tableau_poids[i-1, j]
         if (modele =="logistique_1") :
            tableau_poids[i, j] = x + b_logistique1(x, k, N) * finesse_pas + sigma_logistique1(c, x) * np.random.normal(0, sqrt(finesse_pas))
         if (modele =="muller_feuga") :
            tableau_poids[i, j] = x + b_muller(x, k, N) * finesse_pas + sigma_muller(c, x) * np.random.normal(0, sqrt(finesse_pas))
         if (modele =="logistique_2A") :
            tableau_poids[i, j] = x + b_logi2A(x, k, N) * finesse_pas + sigma_logi2A(c, x) * np.random.normal(0, sqrt(finesse_pas))
         if (modele =="logistique_2B") :
            tableau_poids[i, j] = x + b_logi2B(x, k, N) * finesse_pas + sigma_logi2B(c, x) * np.random.normal(0, sqrt(finesse_pas))
    
    return tableau_poids

# METHODE MILSTEIN
def milstein(N, temps_final, finesse_pas, c, mu, modele):

    nombre_observations = int(temps_final / finesse_pas) + 1
     # Initialisation du tableau
    tableau_poids = np.zeros((nombre_observations, N))
    X_0 = np.random.uniform(0.5, 2, N)  
    tableau_poids[0, :] = X_0

    # Remplissage du tableau
    for i in range(1, nombre_observations):
    
        # On ordonne les poids
        poids_ordonnes = tableau_poids[i-1, :]
        poids_ordonnes = np.sort(poids_ordonnes)

        for j in range(N):
            #k = np.where(poids_ordonnes == tableau_poids[i-1, j])[0][0]
            k = np.searchsorted(poids_ordonnes, tableau_poids[i-1, j])
            delta_Bt = np.random.normal(0, sqrt(finesse_pas))
            x = tableau_poids[i-1, j]
            if (modele == 'logistique_1'):
             tableau_poids[i, j] = x + b_logistique1(x, k, N) * finesse_pas + sigma_logistique1(c, x) * delta_Bt + (1/2)*sigma_logistique1(c, x)*sigma_der_logi1(c,x)*(delta_Bt**2-finesse_pas)
            if (modele == 'muller_feuga'):
             tableau_poids[i, j] = x + b_muller(x, k, N) * finesse_pas + sigma_muller(c, x) * delta_Bt + (1/2)*sigma_muller(c, x)*sigma_der_muller(c,x)*(delta_Bt**2-finesse_pas)
            if (modele == 'logistique_2A'):
             tableau_poids[i, j] = x + b_logi2A(x, k, N) * finesse_pas + sigma_logi2A(c, x) * delta_Bt + (1/2)*sigma_logi2A(c, x)*sigma_der_logi2A(c,x)*(delta_Bt**2-finesse_pas)
            if (modele == 'logistique_2B'):
             tableau_poids[i, j] = x + b_logi2B(x, k, N) * finesse_pas + sigma_logi2B(c, x) * delta_Bt + (1/2)*sigma_logi2B(c, x)*sigma_der_logi2B(c,x)*(delta_Bt**2-finesse_pas)
    return tableau_poids

# METHODE RUNGE KUTTA ORDRE 1.5
def RK(N, temps_final, finesse_pas, c, mu, modele):
    nombre_observations =int(temps_final / finesse_pas) + 1
    # Initialisation du tableau
    tableau_poids = np.zeros((nombre_observations, N))
    X_0 = np.random.uniform(0.5, 2, N)  
    tableau_poids[0, :] = X_0

    # Remplissage du tableau
    for i in range(1, nombre_observations):
    
        # On ordonne les poids
     poids_ordonnes = tableau_poids[i-1, :]
     poids_ordonnes = np.sort(poids_ordonnes)

     for j in range(N):
         k = np.searchsorted(poids_ordonnes, tableau_poids[i-1, j])
         delta_Bt = np.random.normal(0, sqrt(finesse_pas))
         x = tableau_poids[i-1, j]

         if (modele =="logistique_1") :
            x_prime = x + b_logistique1(x,k,N)*finesse_pas + sigma_logistique1(c,x)*sqrt(finesse_pas)
            tableau_poids[i, j] = x + b_logistique1(x, k, N) * finesse_pas + sigma_logistique1(c, x) * delta_Bt + (1/2)*(sigma_logistique1(c, x_prime)-sigma_logistique1(c,x))*(delta_Bt**2-finesse_pas)/sqrt(finesse_pas)
         if (modele =="muller_feuga") :
            x_prime = x + b_muller(x,k,N)*finesse_pas + sigma_muller(c,x)*sqrt(finesse_pas)
            tableau_poids[i, j] = x + b_muller(x, k, N) * finesse_pas + sigma_muller(c, x) * delta_Bt + (1/2)*(sigma_muller(c, x_prime)-sigma_muller(c,x))*(delta_Bt**2-finesse_pas)/sqrt(finesse_pas)
         if (modele =="logistique_2A") :
            x_prime = x + b_logi2A(x,k,N)*finesse_pas + sigma_logi2A(c,x)*sqrt(finesse_pas)
            tableau_poids[i, j] = x + b_logi2A(x, k, N) * finesse_pas + sigma_logi2A(c, x) * delta_Bt + (1/2)*(sigma_logi2A(c, x_prime)-sigma_logi2A(c,x))*(delta_Bt**2-finesse_pas)/sqrt(finesse_pas)
         if (modele =="logistique_2B") :
            x_prime = x + b_logi2B(x,k,N)*finesse_pas + sigma_logi2B(c,x)*sqrt(finesse_pas)
            tableau_poids[i, j] = x + b_logi2B(x, k, N) * finesse_pas + sigma_logi2B(c, x) * delta_Bt + (1/2)*(sigma_logi2B(c, x_prime)-sigma_logi2B(c,x))*(delta_Bt**2-finesse_pas)/sqrt(finesse_pas)
    
    return tableau_poids

