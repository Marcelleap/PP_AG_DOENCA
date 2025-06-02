/* Esse programa é o codigo sequencia relicionado ao trabalho de pp dessa maneira ele contem esplicações de respectivas funções a serem ultilizadas */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define BITS_POR_GENE 8
#define NUM_GENES 5
#define TAM_CROMOSSOMO (BITS_POR_GENE * NUM_GENES)

#define TAM_POPULACAO 50
#define MAX_GER 100

#define PROB_CROSSOVER 0.8
#define PROB_MUTACAO 0.01

// Pesos dos genes para fitness
#define PESO_RAD 0.4
#define PESO_FOME 0.2
#define PESO_ENERGIA 0.4
#define PESO_SANIDADE 0.2
#define PESO_DOENCA 0.2

typedef struct {
    int genes[TAM_CROMOSSOMO];
    float fitness;
} Individuo;

// Função para converter binário em inteiro
int binario_para_inteiro(int *genes, int inicio, int fim) {
    int valor = 0;
    for (int i = inicio; i <= fim; i++) {
        valor = (valor << 1) | genes[i];
    }
    return valor;
}

// Cálculo de fitness
float calcular_fitness(Individuo individuo) {
    int rad = binario_para_inteiro(individuo.genes, 0, 7);
    int fome = binario_para_inteiro(individuo.genes, 8, 15);
    int energia = binario_para_inteiro(individuo.genes, 16, 23);
    int sanidade = binario_para_inteiro(individuo.genes, 24, 31);
    int doenca = binario_para_inteiro(individuo.genes, 32, 39);

    float fitness = (rad * PESO_RAD) + (fome * PESO_FOME) + (energia * PESO_ENERGIA)
                  + (sanidade * PESO_SANIDADE) + (doenca * PESO_DOENCA);

    return fitness;
}

// Inicializar população
void inicializar_populacao(Individuo *populacao) {
    for (int i = 0; i < TAM_POPULACAO; i++) {
        for (int j = 0; j < TAM_CROMOSSOMO; j++) {
            populacao[i].genes[j] = rand() % 2;
        }
        populacao[i].fitness = calcular_fitness(populacao[i]);
    }
}

// Seleção por torneio
Individuo torneio(Individuo *populacao) {
    int a = rand() % TAM_POPULACAO;
    int b = rand() % TAM_POPULACAO;
    if (populacao[a].fitness > populacao[b].fitness) {
        return populacao[a];
    } else {
        return populacao[b];
    }
}

// Cruzamento de um ponto
void cruzamento(Individuo pai1, Individuo pai2, Individuo *filho1, Individuo *filho2) {
    if ((float)rand() / RAND_MAX < PROB_CROSSOVER) {
        int ponto = rand() % TAM_CROMOSSOMO;
        for (int i = 0; i < ponto; i++) {
            filho1->genes[i] = pai1.genes[i];
            filho2->genes[i] = pai2.genes[i];
        }
        for (int i = ponto; i < TAM_CROMOSSOMO; i++) {
            filho1->genes[i] = pai2.genes[i];
            filho2->genes[i] = pai1.genes[i];
        }
    } else {
        *filho1 = pai1;
        *filho2 = pai2;
    }
    filho1->fitness = calcular_fitness(*filho1);
    filho2->fitness = calcular_fitness(*filho2);
}

// Mutação
void mutacao(Individuo *individuo) {
    for (int i = 0; i < TAM_CROMOSSOMO; i++) {
        if ((float)rand() / RAND_MAX < PROB_MUTACAO) {
            individuo->genes[i] = !individuo->genes[i];
        }
    }
    individuo->fitness = calcular_fitness(*individuo);
}

// Impressão
void imprimir_individuo(Individuo ind) {
    for (int i = 0; i < TAM_CROMOSSOMO; i++) {
        printf("%d", ind.genes[i]);
    }
    int rad = binario_para_inteiro(ind.genes, 0, 7);
    int fome = binario_para_inteiro(ind.genes, 8, 15);
    int energia = binario_para_inteiro(ind.genes, 16, 23);
    int sanidade = binario_para_inteiro(ind.genes, 24, 31);
    int doenca = binario_para_inteiro(ind.genes, 32, 39);

    printf(" | Rad: %d, Fome: %d, Energia: %d, Sanidade: %d, Doença: %d | Fitness: %.2f\n",
           rad, fome, energia, sanidade, doenca, ind.fitness);
}

// ------------------ MAIN ------------------
int main() {
    srand(time(NULL));

    Individuo populacao[TAM_POPULACAO];
    Individuo nova_geracao[TAM_POPULACAO];

    inicializar_populacao(populacao);

    for (int geracao = 0; geracao < MAX_GER; geracao++) {
        for (int i = 0; i < TAM_POPULACAO; i += 2) {
            Individuo pai1 = torneio(populacao);
            Individuo pai2 = torneio(populacao);

            Individuo filho1, filho2;
            cruzamento(pai1, pai2, &filho1, &filho2);

            mutacao(&filho1);
            mutacao(&filho2);

            nova_geracao[i] = filho1;
            nova_geracao[i + 1] = filho2;
        }

        // Substituir antiga população
        for (int i = 0; i < TAM_POPULACAO; i++) {
            populacao[i] = nova_geracao[i];
        }

        // Encontrar o melhor da geração
        Individuo melhor = populacao[0];
        for (int i = 1; i < TAM_POPULACAO; i++) {
            if (populacao[i].fitness > melhor.fitness) {
                melhor = populacao[i];
            }
        }

        printf("Geração %d | Melhor Fitness: %.2f\n", geracao + 1, melhor.fitness);
        imprimir_individuo(melhor);
    }

    return 0;
}




