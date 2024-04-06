void pivoter(vector<vector<float>> &A, vector<float> &B, int dim, int pi)
{ // fonction pivotant la matrice A ayant comme pivot la ligne pi
    for (int i = pi + 1; i < dim; i++)
    {
        for (int j = pi + 1; j < dim; j++)
        {
            //A[i][j] = A[i][j] - ( (A[i][pi] / A[pi][pi]) * A[pi][j] )
            A.at(i).at(j) = A.at(i).at(j) - ((A.at(i).at(pi) / A.at(pi).at(pi)) * A.at(pi).at(j));
        }
        B[i] = B[i] - ((A.at(i).at(pi) / A.at(pi).at(pi)) * B.at(pi));
        A.at(i).at(pi) = 0;
    }
    cout << "Apres pivotage, " << endl; // afficher le système d'équations après pivotage
    afficher(A, B, dim);
}