import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from fonctions import euler_explicite, milstein, RK  # et toutes les fonctions qu’ils appellent

def demander_entree():
    N = int(input("Nombre de poissons : "))
    temps_final = float(input("Temps final : "))
    finesse_pas = float(input("Finesse du pas : "))
    c = float(input("Paramètre c : "))
    mu = float(input("Paramètre mu (non utilisé ici, mais conservé) : "))
    modele = input("Modèle (logistique_1 / muller_feuga / logistique_2A / logistique_2B) : ")
    methode = input("Méthode (euler_explicite / milstein / RK): ")
    nb_simulations = int(input("Nombre de simulations : "))
    return N, temps_final, finesse_pas, c, mu, modele, methode, nb_simulations

def tracer_stats(simulations, temps, N):
    poids_min = np.min(simulations, axis=0)
    poids_max = np.max(simulations, axis=0)
    poids_moyen = np.mean(simulations, axis=0)

    # Index du poisson le plus maigre et le plus gros à t=0
    indices_ordre = np.argsort(simulations[0, 0, :])
    idx_maigre = indices_ordre[0]
    idx_gros = indices_ordre[-1]

    plt.figure(figsize=(12, 6))
    plt.fill_between(temps, poids_min[:, idx_maigre], poids_max[:, idx_maigre], color='pink', alpha=0.3, label='Min/Max maigre')
    plt.plot(temps, poids_moyen[:, idx_maigre], color='deeppink', label='Moyenne maigre')

    plt.fill_between(temps, poids_min[:, idx_gros], poids_max[:, idx_gros], color='skyblue', alpha=0.3, label='Min/Max gros')
    plt.plot(temps, poids_moyen[:, idx_gros], color='blue', label='Moyenne gros')

    #plt.fill_between(temps, np.min(poids_min, axis=1), np.max(poids_max, axis=1), color='lightgreen', alpha=0.3, label='Min/Max global')
    plt.plot(temps, np.mean(poids_moyen, axis=1), color='green', label='Moyenne globale')

    plt.xlabel("Temps")
    plt.ylabel("Poids")
    plt.title("Évolution statistique du poids des poissons")
    plt.legend()
    plt.grid(True)
    plt.show()


# Dictionnaire des méthodes disponibles
methodes_disponibles = {
    "euler_explicite": euler_explicite,
    "milstein": milstein,
    "RK": RK
}

def main():
    N, temps_final, finesse_pas, c, mu, modele, methode, nb_simulations = demander_entree()
    nombre_observations = int(temps_final / finesse_pas) + 1
    temps = np.linspace(0, temps_final, nombre_observations)

    simulations = np.zeros((nb_simulations, nombre_observations, N))

    for s in range(nb_simulations):
        if methode in methodes_disponibles:
            tableau = methodes_disponibles[methode](N, temps_final, finesse_pas, c, mu, modele)
        else:
            raise ValueError("Méthode non reconnue. Choisir entre euler_explicite, milstein, RK.")
        simulations[s] = tableau

    tracer_stats(simulations, temps, N)

if __name__ == "__main__":
    main()