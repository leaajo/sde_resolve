Le 17 avril 2025 

-------------------------
🐟 Simulation de Poissons

Permet de simuler l'évolution du poids de poissons selon différents modèles biologiques et méthodes numériques.

Interface graphique réalisée avec Tkinter, calculs avec NumPy et visualisation avec Matplotlib.


--------------------------
📦 INSTALLATION
--------------------------
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt

--------------------------
🚀 LANCER L’APPLICATION
--------------------------

Toujours depuis le terminal, avec l’environnement virtuel activé :

   source venv/bin/activate
   python chemin_vers_le_dossier/new_interface.py

(Remplacer "chemin_vers_le_dossier" par le chemin réel vers le fichier new_interface.py.)

--------------------------
🔧 PARAMÈTRES DE SIMULATION
--------------------------

Dans l’interface, tu peux renseigner les paramètres suivants :

- Nombre de poissons
- Temps final
- Finesse du pas
- Paramètre c
- Paramètre mu
- Modèle biologique (logistique_1, muller_feuga, etc.)
- Méthode numérique (euler_explicite, milstein, RK)
- Nombre de simulations

--------------------------
📊 RÉSULTATS
--------------------------

Un graphique s’affiche à la fin de la simulation, montrant :
- Les courbes min/max/moyenne des poissons les plus maigres et les plus gros
- La moyenne globale du poids au fil du temps


--------------------------
🐍 DÉPENDANCES PRINCIPALES
--------------------------

- numpy
- matplotlib
- tkinter (déjà inclus avec Python standard)
