import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from fonctions import euler_explicite, milstein, RK

# dico des méthodes disponibles
methodes_disponibles = {
    "euler_explicite": euler_explicite,
    "milstein": milstein,
    "RK": RK
}

# Fonction de simulation
def lancer_simulation():
    print("Bouton cliqué")

    try:
        N = int(entry_N.get())
        temps_final = float(entry_temps_final.get())
        finesse_pas = float(entry_finesse_pas.get())
        c = float(entry_c.get())
        mu = float(entry_mu.get())
        modele = combo_modele.get()
        methode = combo_methode.get()
        nb_simulations = int(entry_nb_simulations.get())

        nombre_observations = int(temps_final / finesse_pas) + 1
        temps = np.linspace(0, temps_final, nombre_observations)
        simulations = np.zeros((nb_simulations, nombre_observations, N))

        for s in range(nb_simulations):
            if methode in methodes_disponibles:
                tableau = methodes_disponibles[methode](N, temps_final, finesse_pas, c, mu, modele)
            else:
                raise ValueError("Méthode non reconnue.")
            simulations[s] = tableau

        # Appel de tracer_stats
        tracer_stats(simulations, temps, N)

    except Exception as e:
        messagebox.showerror("Erreur", f"Une erreur est survenue :\n{e}")

# Fonction tracer_stats
def tracer_stats(simulations, temps, N):
    poids_min = np.min(simulations, axis=0)
    poids_max = np.max(simulations, axis=0)
    poids_moyen = np.mean(simulations, axis=0)

    indices_ordre = np.argsort(simulations[0, 0, :])
    idx_maigre = indices_ordre[0]
    idx_gros = indices_ordre[-1]

    plt.figure(figsize=(12, 6))
    plt.fill_between(temps, poids_min[:, idx_maigre], poids_max[:, idx_maigre], color='pink', alpha=0.3, label='Min/Max maigre')
    plt.plot(temps, poids_moyen[:, idx_maigre], color='deeppink', label='Moyenne maigre')

    plt.fill_between(temps, poids_min[:, idx_gros], poids_max[:, idx_gros], color='skyblue', alpha=0.3, label='Min/Max gros')
    plt.plot(temps, poids_moyen[:, idx_gros], color='blue', label='Moyenne gros')

    plt.plot(temps, np.mean(poids_moyen, axis=1), color='green', label='Moyenne globale')

    plt.xlabel("Temps")
    plt.ylabel("Poids")
    plt.title("Évolution statistique du poids des poissons")

    # Texte des paramètres
    param_text = (
        f"Nombre de poissons (N) : {N}\n"
        f"Temps final : {entry_temps_final.get()}\n"
        f"Finesse du pas : {entry_finesse_pas.get()}\n"
        f"Paramètre c : {entry_c.get()}\n"
        f"Paramètre mu : {entry_mu.get()}\n"
        f"Modèle : {combo_modele.get()}\n"
        f"Méthode : {combo_methode.get()}\n"
        f"Nb simulations : {entry_nb_simulations.get()}"
    )

    # Ajouter le texte dans un encadré
    plt.gcf().text(0.72, 0.5, param_text, fontsize=9, bbox=dict(facecolor='white', edgecolor='black'))

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# === Interface ===
fenetre = tk.Tk()
fenetre.title("Simulation de poissons 🐟")
fenetre.geometry("400x500")

# Champs d'entrée
def ajouter_champ(label, row):
    tk.Label(fenetre, text=label).grid(row=row, column=0, sticky="w", padx=10, pady=5)
    champ = tk.Entry(fenetre)
    champ.grid(row=row, column=1)
    return champ

entry_N = ajouter_champ("Nombre de poissons", 0)
entry_temps_final = ajouter_champ("Temps final", 1)
entry_finesse_pas = ajouter_champ("Finesse du pas", 2)
entry_c = ajouter_champ("Paramètre c", 3)
entry_mu = ajouter_champ("Paramètre mu", 4)
entry_nb_simulations = ajouter_champ("Nb simulations", 5)

tk.Label(fenetre, text="Modèle").grid(row=6, column=0, sticky="w", padx=10, pady=5)
combo_modele = ttk.Combobox(fenetre, values=["logistique_1", "muller_feuga", "logistique_2A", "logistique_2B"])
combo_modele.current(0)
combo_modele.grid(row=6, column=1)

tk.Label(fenetre, text="Méthode").grid(row=7, column=0, sticky="w", padx=10, pady=5)
combo_methode = ttk.Combobox(fenetre, values=["euler_explicite", "milstein", "RK"])
combo_methode.current(0)
combo_methode.grid(row=7, column=1)

# Bouton lancer
btn_simuler = tk.Button(fenetre, text="Lancer la simulation", command=lancer_simulation, bg="blue")
btn_simuler.grid(row=8, columnspan=2, pady=20)

fenetre.mainloop()
