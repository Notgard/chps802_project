#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

void afficher(vector<vector<float>> A, vector<float> B, int dim);
float max(float a, float b);
int plusGrandPivot(vector<vector<float>> A, int dim, int pi);
void pivoter(vector<vector<float>> &A, vector<float> &B, int dim, int pi);
void permuter(vector<vector<float>> &A, vector<float> &B, int dim, int L1, int L2);
void triangulariser(vector<vector<float>> &A, vector<float> &B, int dim);
void resoudre(vector<vector<float>> A, vector<float> B, int dim, vector<float> X);
void recupererFichier(vector<vector<float>> &A, vector<float> &B, int &dim);
vector<float> entrerVectorLigne(string ligne);

int main()
{
    cout << "*****Méthode d'élimination de Gauss avec pivotage partiel*****" << endl;
    int dim = 4;
    vector<vector<float>> A; // initialisation et récupération des données
    vector<float> B;
    recupererFichier(A, B, dim);
    vector<float> X;
    for (int i = 0; i < dim; i++)
        X.push_back(0);

    cout << "Voici le système sous forme matriciel :" << endl; // les étapes de résolution
    afficher(A, B, dim);
    triangulariser(A, B, dim);
    resoudre(A, B, dim, X);

    return 0;
}

void afficher(vector<vector<float>> A, vector<float> B, int dim)
{ // affichage du système sous forme matriciel
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            cout << A.at(i).at(j) << "\t";
        }
        cout << "\tx" << i + 1 << "\t\t" << B.at(i) << endl;
    }
}
float max(float a, float b)
{ // fonction retournant la valeur du nombre plus grand par rapport à 0
    return fabs(a) < fabs(b) ? b : a;
}
int plusGrandPivot(vector<vector<float>> A, int dim, int pi)
{ // fonction retournant le plus grand pivot dans une colonne
    int Lpiv = 0;
    float piv = 0;
    float a = A.at(pi).at(pi);
    for (int i = pi + 1; i < dim; i++)
    {
        piv = max(a, A.at(i).at(pi)); // trouver le plus grand pivot
        a = piv;
    }
    for (int i = 0; i < dim; i++)
    {
        if (A.at(i).at(pi) == piv)
        { // trouver l'adresse du plus grand pivot
            Lpiv = i;
        }
    }
    cout << endl
         << "Le plus grand pivot est " << piv << " et se trouve à la ligne " << Lpiv + 1 << endl;

    return Lpiv;
}
void pivoter(vector<vector<float>> &A, vector<float> &B, int dim, int pi)
{ // fonction pivotant la matrice A ayant comme pivot la ligne pi
    for (int i = pi + 1; i < dim; i++)
    {
        for (int j = pi + 1; j < dim; j++)
        {
            A.at(i).at(j) = A.at(i).at(j) - ((A.at(i).at(pi) / A.at(pi).at(pi)) * A.at(pi).at(j));
        }
        B[i] = B[i] - ((A.at(i).at(pi) / A.at(pi).at(pi)) * B.at(pi));
        A.at(i).at(pi) = 0;
    }
    cout << "Apres pivotage, " << endl; // afficher le système d'équations après pivotage
    afficher(A, B, dim);
}
void permuter(vector<vector<float>> &A, vector<float> &B, int dim, int L1, int L2)
{ // fonction permutant les lignes L1 et L2 de chacun des matrices A et B
    cout << "Permutation de la ligne " << L1 + 1 << " avec la ligne " << L2 + 1 << " : " << endl;
    for (int i = 0; i < dim; i++)
    {
        float t = A.at(L1).at(i); // permuter les lignes de la matrice A
        A.at(L1).at(i) = A.at(L2).at(i);
        A.at(L2).at(i) = t;
    }
    float t2 = B.at(L1); // puis permuter les lignes de la matrice B
    B.at(L1) = B.at(L2);
    B.at(L2) = t2;
}
void triangulariser(vector<vector<float>> &A, vector<float> &B, int dim)
{ // fonction transformant la matrice A sous la forme triangulaire superieure
    for (int i = 0; i < (dim - 1); i++)
    {
        // etape 1 : initialisation
        int pi = i;

        // etape 2 : recherche du plus grand pivot
        int Lpiv = plusGrandPivot(A, dim, pi);

        // etape 3 : ramener la ligne où se trouve le plus grand pivot trouve en position de pivot
        if (pi != Lpiv)
        { // permuter seulement si ce ne sont pas les memes lignes
            cout << endl;
            permuter(A, B, dim, pi, Lpiv);
            afficher(A, B, dim);
            cout << endl;
        }

        // etape 4 : elimination
        pivoter(A, B, dim, pi);
    }
}
void resoudre(vector<vector<float>> A, vector<float> B, int dim, vector<float> X)
{ // fonction résolvant la forme finale du systeme
    int j = 0;
    float result = 0;
    for (int i = dim - 1; i >= 0; i--)
    { // commencer à la dernière ligne
        for (result = 0, j = i; j < dim; j++)
        {
            cout << "j: " << j << ", dim: " << dim << "i:" << i << endl;
            if (i != j) {
                result += (A.at(i).at(j)) * X.at(j);
                cout << "in:" << A.at(i).at(j) << " " << X.at(j) << endl;
            }
            cout << "what: " << result << endl;
        }
        cout << i << " X[i]=" <<  B.at(i) 
        X.at(i) = (B.at(i) - result) / A.at(i).at(i);
    }
    cout << endl
         << "La solution est :" << endl;
    for (int i = 0; i < dim; i++)
    { // afficher la solution
        cout << "x" << i + 1 << "=" << X.at(i) << endl;
    }
}
void recupererFichier(vector<vector<float>> &A, vector<float> &B, int &dim)
{ // fonction récupérant les données dans le fichier
    fstream fichier;
    string ligne = "";
    int nombreLigne = 0;
    fichier.open("testcpp.txt", ios::in); // ouvrir le fichier "data.txt"
    if (!fichier)
    {
        cout << "Erreur lors de lecture de fichier" << endl; // si le fichier n'existe pas ou ne se trouve pas dans le dossier actuel
        cout << "Le fichier n'existe pas ou ne se trouve pas dans le dossier actuel" << endl;
        exit(-1); /// afficher erreur et arrêter le programme
    }
    while (getline(fichier, ligne))
    {
        if (nombreLigne == 0)
            dim = stoi(ligne); // affecter la valeur de la ligne 0 à la variable dim
        else if (nombreLigne < dim + 1)
            A.push_back(entrerVectorLigne(ligne)); // affecter les lignes suivantes de longueur dim à la matrice A
        else
        {
            B.push_back(stof(ligne)); // affecter le reste à la matrice B
        }
        nombreLigne++;
    }
    fichier.close();
}
vector<float> entrerVectorLigne(string ligne)
{ // fonction récupérant les valeurs dans une ligne
    int i = 0;
    vector<float> retour;
    ligne += " "; // ajouter un string " " à la ligne pour avoir tous les string de la ligne
    int longueur = ligne.length();
    string chiffreString = "";
    while (i < longueur)
    {
        if ((ligne.at(i) >= '0' && ligne.at(i) <= '9') || ligne.at(i) == '-' || ligne.at(i) == '.')
            chiffreString += ligne.at(i); // récupérer seulement les chiffres et les signes
        else
        {
            if (chiffreString.length() > 0)
                retour.push_back(stof(chiffreString)); // si ce n'est pas un chiffre ou un signe, affecter la valeur au tableau retour
            chiffreString = "";                        /// puis vider la variable
        }
        i++;
    }
    return retour;
}