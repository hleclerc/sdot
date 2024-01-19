* trouver les cellules associées aux points et au edges infinis
* faire la version robuste


La proposition pour la création de bord à partir des cellules infinies, c'est de considérer que c'est un convex hull.

Pour la direction, si on peut trouver un point strictement à l'intérieur du domaine, on l'utilise ce point comme départ
  Si on n'est pas capables de trouver un point, 

De base, chaque cellule infinie donne un point dans les bords de la transformée de Legendre
  C'est possible d'avoir 3 points alignés
  Rq: on pourrait faire une matrice de covariance pour voir les directions "bloquées", mais ça va poser des problèmes de précision numérique
  On pourrait aussi envisager de faire des calculs exacts
  Dans ce cas, on pourrait trouver les vecteurs propres de la matrice de covariance
    => les vps nulles vont donner des égalitées
    => les vps non nulles vont donner des inégalitées, et on sera capables de donner le sens avec un point intérieur

Pb: avec 2 cellules, les frontières se "replient".
  On pourrait imaginer travailler dans un espace de plus petite dimension lorsque le noyau de la matrice de covarience est non vide


Pourquoi faudrait-il faire les opérations en précision infinie ?
  On aimerait bien repérer les cas "aplatis" où les gradients et les bords se situent dans un sous-domaine
  On pourrait commencer à 

Rq: les cas aplatis arrivent lorsque le span des coupes n'occupe qu'une partie de l'espace
  Prop: dans un premier temps, on fait simplement les calculs avec des flottants. Ça marchera lorsque le nombre de point n'est pas suffisant.

On cherche un point intérieur dans la transformée de Legendre pour vérifier que boundaries sont bien orientés.
  Les bords viennent des cellules infinies qui donnent à chaque fois un points dans la transformée
  S'il n'y avait que de fonctions affines (pas de bords), il suffirait de faire la moyenne des directions
  Rq: si on a 2 bords en 1D, il n'y a pas de cellule infinie
  Si on a 1 bord et une cellule, ...

Prop: on cherche le point intérieur en travaillant sur les vertices de diagramme initial. En 1D
  * si on a 2 fonctions affines, on prend la direction moyenne
  * si on a 1 fonction affine et un bord, on prend 

# en 1D
max(-2.0*y_0 - 5.0, 0, 3.0*y_0 - 4.0)
max(-2.5*y_0, 1.33333333333333*y_0) for -1.0*y_0 <= 2.0, 1.0*y_0 <= 3.0

# en 2D
max(-2.5*y_0, 1.33333333333333*y_0) for -1.0*y_0 <= 2.0, -1.0*y_1 <= 0.0, 1.0*y_0 <= 3.0, 1.0*y_1 <= 0.0
