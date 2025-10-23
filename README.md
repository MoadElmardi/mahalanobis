# Description

Ce projet est une implémentation du calcul de la distance de Mahalanobis entre un vecteur et une matrice. Le type de nombre utilisé est les entiers signés. 

Le calcul est bien implémenté, mais les résultats ne sont pas satisfaisants. En effet, le calcul de la distance de Mahalanobis contient beaucoup de divisions. Les résultats de ces divisions sont généralement des flottants. Notre utilisation des entiers chiffrés diminue la précision avec chaque division. Ainsi, à la fin du calcul, la distance de Mahalanobis est généralement égale à 0.

Une solution pour ce problème est d'utiliser la librairie des flottants implémentée dans le projet _librairie-flottants_ pour faire le calcul avec le type des flottans au lieu des entiers signés.