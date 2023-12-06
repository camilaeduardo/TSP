#include "Data.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <time.h>
#include <ctime>

using namespace std;

typedef struct Solution{
    vector<int> sequence;
    double cost;
} Solution;

typedef struct insertionInfo{
    int noInserido;
    int arestaRemovida;
    double custo;
} insertionInfo;

void showSolution(Solution *s) {
    for(int i = 0; i < s->sequence.size() - 1; i ++)
        cout << s->sequence[i] << "-->";
    cout << s->sequence.back() << endl;
}

void calculaCusto(Solution *s, Data &data) {
    s->cost = 0;

    for(int i = 0; i < s->sequence.size() - 1; i++) {
        s->cost += data.getDistance(s->sequence[i], s->sequence[i+1]);
    }
}

vector<int> escolhe3NosAleatorios(vector<int> *cl, int dimension) {
    vector<int> v, sequence;

    for(int i = 2; i <= dimension; i++)
        v.push_back(i);

    sequence.push_back(1);
    sequence.push_back(1);

    for(int  i = 0; i < 3; i++) {
        int index = rand() % v.size();
        sequence.insert(sequence.end() - 1,v[index]);
        v.erase(v.begin() + index);
    }

    *cl = v;
    return sequence;
}

vector <insertionInfo> calcularCustoInsercao(Solution *s, vector<int> *cl, Data &data) {
    vector<insertionInfo> custoInsercao((s->sequence.size() - 1) * cl->size());

    int l = 0;

    for (int a = 0; a < s->sequence.size() - 1; a++) {
        int i = s->sequence[a];
        int j = s->sequence[a + 1];
        for(int x = 0; x < cl->size(); x++) {
            custoInsercao[l].custo = data.getDistance(i, (*cl)[x]) + data.getDistance((*cl)[x], j) - data.getDistance(i, j);
            custoInsercao[l].noInserido = x;
            custoInsercao[l].arestaRemovida = a;
            l++;
        }
    }

    return custoInsercao;
}

bool verificaOrdem(insertionInfo x1, insertionInfo x2) {
    return (x1.custo < x2.custo);
}

vector <insertionInfo> ordernarEmOrdemCrescente(vector<insertionInfo> custoInsercao) {
    sort(custoInsercao.begin(), custoInsercao.end(), verificaOrdem);
    return custoInsercao;
}

void inserirNaSolucao(Solution *s, insertionInfo custoInsercao, vector<int> *cl) {
   s->sequence.insert(s->sequence.begin() + custoInsercao.arestaRemovida + 1, (*cl)[custoInsercao.noInserido]);

   cl->erase(cl->begin() + custoInsercao.noInserido);
}

Solution Construcao(Data &data) {
    Solution s;
    vector<int> cl;
    s.sequence = escolhe3NosAleatorios(&cl, data.getDimension());

    while(!cl.empty()) {
        vector<insertionInfo> custoInsercao = calcularCustoInsercao(&s, &cl, data);
        custoInsercao = ordernarEmOrdemCrescente(custoInsercao);
        double alpha = (double) rand() / RAND_MAX;
        int selecionado = rand() % ((int) ceil(alpha * custoInsercao.size()));
        inserirNaSolucao(&s, custoInsercao[selecionado], &cl);

    }

    return s;
}

void swap(Solution *s, int best_i, int best_j) {
    int value_i = s->sequence[best_i];
    s->sequence[best_i] = s->sequence[best_j];
    s->sequence[best_j] = value_i;
}

bool bestImprovementSwap(Solution *s, Data &data) {
    double best_delta = 0;
    int best_i, best_j;

    for(int i = 1; i < s->sequence.size() - 1; i++) {
        int vi = s->sequence[i];
        int vi_next = s->sequence[i+1];
        int vi_prev = s->sequence[i-1];
        for(int j = 1; j < s->sequence.size() - 1; j++) {
            if(i == j) continue;
            int vj = s->sequence[j];
            int vj_next = s->sequence[j+1];
            int vj_prev = s->sequence[j-1];
            double delta;
            if(j == (i + 1)) {
                delta = data.getDistance(vi_prev, vj) + data.getDistance(vj, vi) + data.getDistance(vi, vj_next)
                - data.getDistance(vi_prev, vi) - data.getDistance(vi, vj) - data.getDistance(vj, vj_next);
            }
            else if(j == (i - 1)) {
                    delta = data.getDistance(vj_prev, vi) + data.getDistance(vi, vj) + data.getDistance(vj, vi_next)
                    - data.getDistance(vj_prev, vj) - data.getDistance(vj, vi) - data.getDistance(vi, vi_next);
            }
                else {
                    delta = data.getDistance(vi_prev, vj) + data.getDistance(vj, vi_next) + data.getDistance(vj_prev, vi) + data.getDistance(vi, vj_next)
                    - data.getDistance(vi_prev, vi) - data.getDistance(vi, vi_next) - data.getDistance(vj_prev, vj) - data.getDistance(vj, vj_next);

                }
            if(delta < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }
    if(best_delta < 0) {
        s->cost = s->cost + best_delta;
        swap(s, best_i, best_j);
        return true;
    }
    return false;
}

void reinsertion(Solution *s, int best_i, int best_j) {
    int value_j = s->sequence[best_j];

    if(best_i > best_j) {
        s->sequence.insert(s->sequence.begin() + best_i + 1, value_j);
        s->sequence.erase(s->sequence.begin() + best_j);
    }
    else{
        s->sequence.erase(s->sequence.begin() + best_j);
        s->sequence.insert(s->sequence.begin() + best_i + 1, value_j);
    }
}

bool bestImprovementReinsertion(Solution *s, Data &data) {
    double best_delta = 0;
    int best_i, best_j;
    for(int i = 0; i < s->sequence.size() - 1; i++) {
        int vi_1 = s->sequence[i];
        int vi_2 = s->sequence[i + 1];
        for(int j = 1; j < s->sequence.size() - 1; j++) {
            if(j == i || j == (i+1)) continue;
            int vj = s->sequence[j];
            int vj_prev = s->sequence[j - 1];
            int vj_next = s->sequence[j + 1];

            double delta = data.getDistance(vj_prev, vj_next) + data.getDistance(vi_1, vj) + data.getDistance(vj, vi_2)
            - data.getDistance(vj_prev, vj) - data.getDistance(vj, vj_next) - data.getDistance(vi_1, vi_2);

            if(delta < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }
    if(best_delta < 0) {
        s->cost += best_delta;
        reinsertion(s, best_i, best_j);
        return true;
    }

    return false;
}

void orOpt2(Solution *s, int best_i, int best_j) {
    int value1_i = s->sequence[best_i];
    int value2_i = s->sequence[best_i + 1];

    s->sequence.erase(s->sequence.begin() + best_i);
    s->sequence.erase(s->sequence.begin() + best_i);

    if(best_i > best_j) {
        s->sequence.insert(s->sequence.begin() + best_j, value1_i);
        s->sequence.insert(s->sequence.begin() + best_j + 1, value2_i);
    }
    else {
        s->sequence.insert(s->sequence.begin() + best_j - 2, value1_i);
        s->sequence.insert(s->sequence.begin() + best_j - 1, value2_i);
    }
}

bool bestImprovementOrOpt2(Solution *s, Data &data) {
    double best_delta = 0;
    int best_i, best_j;

    for(int i = 1; i < s->sequence.size() - 2; i++) {
        int vi_1 = s->sequence[i];
        int vi_2 = s->sequence[i + 1];
        int vi_prev = s->sequence[i - 1];
        int vi_next = s->sequence[i + 2];

        for(int j = 1; j < s->sequence.size(); j++) {
            if(j == i || j == (i+1) || j ==(i+2)) continue;
            int vj_prev = s->sequence[j - 1];
            int vj_next = s->sequence[j];

            double delta = data.getDistance(vi_prev, vi_next) + data.getDistance(vj_prev, vi_1) + data.getDistance(vi_2, vj_next)
            - data.getDistance(vi_prev, vi_1) - data.getDistance(vi_2, vi_next) - data.getDistance(vj_prev, vj_next);

            if(delta < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if(best_delta < 0) {                                                                                         
        s->cost += best_delta;
        orOpt2(s, best_i, best_j);
        return true;
    }

    return false;
}

void orOpt3(Solution *s, int best_i, int best_j) {
    int value1_i = s->sequence[best_i];
    int value2_i = s->sequence[best_i + 1];
    int value3_i = s->sequence[best_i +2];
    for(int i = 0; i <= 2; i++) {
        s->sequence.erase(s->sequence.begin() + best_i);
    }

    if(best_i > best_j) {
        s->sequence.insert(s->sequence.begin() + best_j, value1_i);
        s->sequence.insert(s->sequence.begin() + best_j + 1, value2_i);
        s->sequence.insert(s->sequence.begin() + best_j + 2, value3_i);
    }
    else{
        s->sequence.insert(s->sequence.begin() + best_j - 3, value1_i);
        s->sequence.insert(s->sequence.begin() + best_j - 2, value2_i);
        s->sequence.insert(s->sequence.begin() + best_j - 1, value3_i);
    }
}

bool bestImprovementOrOpt3(Solution *s, Data &data) {
    double best_delta = 0;
    int best_i, best_j;

    for(int i = 1; i < s->sequence.size() - 3; i++) {
        int vi_1 = s->sequence[i];
        int vi_3 = s->sequence[i + 2];
        int vi_prev = s->sequence[i - 1];
        int vi_next = s->sequence[i + 3];

        for(int j = 1; j < s->sequence.size(); j++) {
            if(j == i || j == i + 1 || j == i + 2 || j == i + 3) continue;
            int vj_prev = s->sequence[j - 1];
            int vj_next = s->sequence[j];

            double delta = data.getDistance(vi_prev, vi_next) + data.getDistance(vj_prev, vi_1) + data.getDistance(vi_3, vj_next)
            - data.getDistance(vi_prev, vi_1) - data.getDistance(vi_3, vi_next) - data.getDistance(vj_prev, vj_next);

            if(delta < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }
    if(best_delta < 0) {
        s->cost += best_delta;
        orOpt3(s, best_i, best_j);
        return true;
    }
    return false;
}

void twoOpt(Solution *s, int best_i, int best_j) {
    int value1_j = s->sequence[best_j];
    int value2_j = s->sequence[best_j + 1];

    s->sequence[best_j] = s->sequence[best_i + 1];
    s->sequence[best_j + 1] = s->sequence[best_i];
    s->sequence[best_i] = value2_j;
    s->sequence[best_i + 1] = value1_j;
}

bool bestImprovement2opt(Solution *s, Data &data) {
    double best_delta = 0;
    int best_i, best_j;

    for(int i = 1; i < s->sequence.size() - 2; i++) {
        int vi_1 = s->sequence[i];
        int vi_2 = s->sequence[i+1];
        int vi_prev = s->sequence[i-1];
        int vi_next = s->sequence[i+2];

        for(int j = 1; j < s->sequence.size() - 2; j++) {
            if(j == i || j == i - 1 || j == i + 1) continue;
            int vj_1 = s->sequence[j];
            int vj_2 = s->sequence[j+1];
            int vj_prev = s->sequence[j-1];
            int vj_next = s->sequence[j+2];

            double delta;

            if(j == i + 2) {
                delta = data.getDistance(vi_prev, vj_2) + data.getDistance(vj_2, vj_1) + data.getDistance(vj_1, vi_2) + data.getDistance(vi_2, vi_1)
                + data.getDistance(vi_1, vj_next) - data.getDistance(vi_prev, vi_1) - data.getDistance(vi_1, vi_2) - data.getDistance(vi_2, vj_1)
                - data.getDistance(vj_1, vj_2) - data.getDistance(vj_2, vj_next);
            }
            else if(j == i - 2) {
                    delta = data.getDistance(vj_prev, vi_2) + data.getDistance(vi_2, vi_1) + data.getDistance(vi_1, vj_2) + data.getDistance(vj_2, vj_1) +
                    data.getDistance(vj_1, vi_next) - data.getDistance(vj_prev, vj_1) - data.getDistance(vj_1, vj_2) - data.getDistance(vj_2, vi_1) -
                    data.getDistance(vi_1, vi_2) - data.getDistance(vi_2, vi_next); 
            }
                else {
                    delta = data.getDistance(vi_prev, vj_2) + data.getDistance(vj_1, vi_next) + data.getDistance(vj_prev, vi_2) + data.getDistance(vi_1, vj_next)
                    + data.getDistance(vi_2, vi_1) + data.getDistance(vj_2, vj_1) - data.getDistance(vi_prev, vi_1) - data.getDistance(vi_1, vi_2) -
                    data.getDistance(vi_2, vi_next) - data.getDistance(vj_prev, vj_1) - data.getDistance(vj_1, vj_2) - data.getDistance(vj_2, vj_next);
                }

            if(delta < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if(best_delta < 0) {
        s->cost += best_delta;
        twoOpt(s, best_i, best_j);
        return true;
    }
    return false;
}

void buscaLocal(Solution *s, Data &data) {
    vector<int> nl = {1,2,3,4,5};
    bool improved = false;
    while(nl.empty() == false) {
        int n = rand() % nl.size();
        switch(nl[n]) {
            case 1: improved = bestImprovementSwap(s, data);
            break;
            case 2: improved = bestImprovement2opt(s, data);
            break;
            case 3: improved = bestImprovementReinsertion(s, data);
            break;
            case 4: improved = bestImprovementOrOpt2(s, data);
            break;
            case 5: improved = bestImprovementOrOpt3(s, data);
            break;
        }

        if(improved)
            nl = {1,2,3,4,5};
        else
            nl.erase(nl.begin() + n);
    }
}

Solution perturbacao(Solution s, Data &data) {
    int tamanho = data.getDimension();
    int size1, size2;
    int tamanho_arr = ceil(tamanho/10.0);

    if(tamanho <= 20) {
        size1 = 2;
        size2 = 2;
    }
    else {
        size1 = (rand() % (tamanho_arr - 1)) + 2;
        size2 = (rand() % (tamanho_arr - 1)) + 2;
    }

    int pos1 = (rand() % (tamanho - size1 - 2)) + 1;
    int pos2;
    int insercao2;

    if(pos1 - size2 >= 1 && pos1 + size1 <= tamanho - size2 - 1)
        insercao2 = (rand() % 2);
    else if(pos1 - size2 >= 1)
        insercao2 = 0;
        else
            insercao2 = 1;

    if(insercao2 == 0)
        pos2 = (rand() % (pos1 - size2)) + 1;
    else
        pos2 = (rand() % (tamanho - size2 - pos1 - size1)) + pos1 + size1;

     int vi_prev = s.sequence[pos1 - 1];
     int vi_next = s.sequence[pos1 + size1];
     int vi_1 = s.sequence[pos1];
     int vi_n = s.sequence[pos1 + size1 - 1];

     int vj_prev = s.sequence[pos2 - 1];
     int vj_next = s.sequence[pos2 + size2];
     int vj_1 = s.sequence[pos2];
     int vj_n = s.sequence[pos2 + size2 - 1];

     double delta;

     if(pos2 == pos1 + size1)
        delta = data.getDistance(vi_prev, vj_1) + data.getDistance(vj_n, vi_1) + data.getDistance(vi_n, vj_next)
        - data.getDistance(vi_prev, vi_1) - data.getDistance(vi_n, vj_1) - data.getDistance(vj_n, vj_next);

     else if(pos1 == pos2 + size2)
            delta = data.getDistance(vj_prev, vi_1) + data.getDistance(vi_n, vj_1) + data.getDistance(vj_n, vi_next)
            - data.getDistance(vj_prev, vj_1) - data.getDistance(vj_n, vi_1) - data.getDistance(vi_n, vi_next);

        else
            delta = data.getDistance(vi_prev, vj_1) + data.getDistance(vj_n, vi_next) + data.getDistance(vj_prev, vi_1) + data.getDistance(vi_n, vj_next)
            - data.getDistance(vi_prev, vi_1) - data.getDistance(vi_n, vi_next) - data.getDistance(vj_prev, vj_1) - data.getDistance(vj_n, vj_next);

    s.cost += delta;

    if(pos1 < pos2) {
        int values_2[size2];

        for(int i = 0; i < size2; i++)
            values_2[i] = s.sequence[pos2 + i];

        for(int i = 0; i < size2; i++)
            s.sequence.erase(s.sequence.begin() + pos2);

        for(int i = 0; i < size1; i++)
            s.sequence.insert(s.sequence.begin() + pos2 + i, s.sequence[pos1 + i]);

        for(int i = 0; i < size1; i++)
            s.sequence.erase(s.sequence.begin() + pos1);

        for(int i = 0; i < size2; i++)
            s.sequence.insert(s.sequence.begin() + pos1 + i, values_2[i]);
    }

    else {
        int values_1[size1];

        for(int i = 0; i < size1; i++)
            values_1[i] = s.sequence[pos1 + i];

        for(int i = 0; i < size1; i++)
            s.sequence.erase(s.sequence.begin() + pos1);

        for(int i = 0; i < size2; i++)
            s.sequence.insert(s.sequence.begin() + pos1 + i, s.sequence[pos2 + i]);

        for(int i = 0; i < size2; i++)
            s.sequence.erase(s.sequence.begin() + pos2);

        for(int i = 0; i < size1; i++)
            s.sequence.insert(s.sequence.begin() + pos2 + i, values_1[i]);
    }

    return s;
}

Solution ILS(int maxIter, int maxIterIls, Data &data) {
    Solution bestOfAll;
    bestOfAll.cost = INFINITY;
    
    for(int i = 0; i < maxIter; i++) {
        Solution s = Construcao(data);
        calculaCusto(&s, data);
        Solution best = s;

        int iterIls = 0;

        while(iterIls <= maxIterIls) {
            buscaLocal(&s, data);
            if(s.cost < best.cost) {
                best = s;
                iterIls = 0;
            }
            s = perturbacao(best, data);
            iterIls++;
        }
        if(best.cost < bestOfAll.cost)
            bestOfAll = best;
    }

    return bestOfAll;
}

int main(int argc, char** argv) {

    srand(time(0));
    auto data = Data(argc, argv[1]);
    data.read();

    int n = data.getDimension();
    int maxIter = 50;
    int maxIterIls;

    if(n >= 150)
        maxIterIls = n/2;
    else
        maxIterIls = n;
    
    clock_t start = clock();

    Solution s = ILS(maxIter, maxIterIls, data);

    clock_t end = clock();
    double time = ((double) (end - start)) / CLOCKS_PER_SEC;
    showSolution(&s);
    cout << "Custo: " << s.cost << endl;
    cout << "Tempo: " << time << endl;

    return 0;
}