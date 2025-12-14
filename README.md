# Générateur de Géométrie de Pale d'Éolienne (Méthode Hansen)

Ce projet contient un script Python (`Palme_nacaXXXX.py`) permettant de dimensionner et de générer la géométrie 3D d'une pale .

Le code s'appuie sur la théorie **BEM (Blade Element Momentum)** et utilise les équations polynomiales de **M.O.L. Hansen** pour calculer la distribution optimale du vrillage (twist) et de la longueur de corde (chord).

## 📋 Fonctionnalités

* **Profil Aérodynamique** : Génération automatique des coordonnées pour les profils NACA 4 chiffres (par défaut **NACA 4418**).
* **Calcul BEM** : Optimisation aérodynamique basée sur le TSR (Tip Speed Ratio) et la finesse maximale.
* **Sécurités** : Prise en compte de cordes minimales et maximales pour la faisabilité de fabrication.
* **Export CAO** : Génération de fichiers de nuages de points (`.txt`) pour chaque section, prêts pour l'importation dans SolidWorks, CATIA, Fusion 360, etc.

## ⚙️ Paramètres du Projet

Les paramètres sont définis en tête du script `Palme_nacaXXXX.py` et peuvent être ajustés :

| Paramètre | Valeur (Défaut) | Description |
| :--- | :--- | :--- |
| `R_rotor` | **1.7 m** | Rayon total (Diamètre = 3.4m) |
| `R_hub` | **0.30 m** | Rayon du moyeu (racine de la pale) |
| `N_blades` | **3** | Nombre de pales |
| `v_inf` | **3 m/s** | Vitesse du vent de conception |
| `w` | **10.5 rad/s** | Vitesse de rotation (~100 RPM) |
| `TSR` | *Calculé (~5.95)* | Tip Speed Ratio ($\lambda$) |
| `NACA_Code` | **"4418"** | Profil utilisé (Courbure/Épaisseur) |
