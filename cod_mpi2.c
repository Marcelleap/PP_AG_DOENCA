#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define BITS_POR_GENE 8
#define NUM_GENES 5
#define TAM_POPULACAO 5
#define TAM_CROMOSSOMO (BITS_POR_GENE * NUM_GENES)
#define MAX_GER 10
#define RAD_MAX 2

#define NUM_ILHAS 2
#define INTERVALO_MIGRACAO 3

float PROB_CROSSOVER;
float PROB_MUTACAO;
int LIMIAR_MORTE;
float PESO_RAD, PESO_FOME, PESO_ENERGIA, PESO_SANIDADE, PESO_DOENCA;

// Struct do indivíduo
typedef struct {
    int genes[TAM_CROMOSSOMO];
    float fitness;
} Individuo;

int binario_para_inteiro(int genes[], int inicio, int fim);
float calcular_fitness(Individuo individuo);
void inicializar_populacao(Individuo populacao[]);
Individuo torneio(Individuo populacao[]);
void cruzamento(Individuo pai, Individuo mae, Individuo *filho1, Individuo *filho2);
Individuo mutacao(Individuo individuo);
void imprimir_populacao(Individuo populacao[], int geracao, int ilha);
void imprimir_individuo(Individuo ind);
void aplicar_morte(Individuo populacao[], int geracao);
void inicializar_parametros(void);

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 3) {
        if (rank == 0) printf("Este programa requer exatamente 3 processos (1 master + 2 ilhas).\n");
        MPI_Finalize();
        return 1;
    }

    double start_time = 0, end_time = 0;
    if (rank == 0) start_time = MPI_Wtime();

    srand(time(NULL) + rank * 1000);
    if (rank != 0) inicializar_parametros();

    Individuo populacao[TAM_POPULACAO];
    Individuo nova_geracao[TAM_POPULACAO];
    Individuo melhor_local, migrante;

    MPI_Datatype MPI_INDIVIDUO;
    int blocklengths[2] = {TAM_CROMOSSOMO, 1};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(Individuo, genes);
    offsets[1] = offsetof(Individuo, fitness);
    MPI_Datatype types[2] = {MPI_INT, MPI_FLOAT};
    MPI_Type_create_struct(2, blocklengths, offsets, types, &MPI_INDIVIDUO);
    MPI_Type_commit(&MPI_INDIVIDUO);

    if (rank == 0) {
        // MASTER: coordena migração
        for (int geracao = 0; geracao < MAX_GER; geracao++) {
            if (geracao % INTERVALO_MIGRACAO == 0) {
                Individuo melhores[NUM_ILHAS];
                // Recebe melhores das ilhas
                for (int ilha = 0; ilha < NUM_ILHAS; ilha++) {
                    MPI_Recv(&melhores[ilha], 1, MPI_INDIVIDUO, ilha+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                // Decide migração circular
                // Ordena índices por fitness decrescente
                int idx[NUM_ILHAS] = {0, 1};
                if (melhores[1].fitness > melhores[0].fitness) {
                    int tmp = idx[0]; idx[0] = idx[1]; idx[1] = tmp;
                }
                Individuo migrantes[NUM_ILHAS];
                migrantes[0] = melhores[idx[1]]; // 2º melhor vai para ilha 0
                migrantes[1] = melhores[idx[0]]; // 1º melhor vai para ilha 1
                // Envia migrantes
                for (int ilha = 0; ilha < NUM_ILHAS; ilha++) {
                    MPI_Send(&migrantes[ilha], 1, MPI_INDIVIDUO, ilha+1, 1, MPI_COMM_WORLD);
                }
                // Log
                printf("\n[MASTER] Migração geração %d:\n", geracao);
                for (int ilha = 0; ilha < NUM_ILHAS; ilha++) {
                    printf("Enviando migrante para ilha %d: fitness %.2f\n", ilha+1, migrantes[ilha].fitness);
                }
            }
        }
        end_time = MPI_Wtime();
        printf("\nTEMPO DE EXECUÇÃO MPI: %.6f segundos\n", end_time - start_time);
    } else {
        // ILHAS: evoluem população
        inicializar_populacao(populacao);
        for (int geracao = 0; geracao < MAX_GER; geracao++) {
            // Evolução local
            for (int i = 0; i < TAM_POPULACAO; i += 2) {
                Individuo pai = torneio(populacao);
                Individuo mae = torneio(populacao);
                Individuo filho1, filho2;
                cruzamento(pai, mae, &filho1, &filho2);
                filho1 = mutacao(filho1);
                filho2 = mutacao(filho2);
                nova_geracao[i] = filho1;
                if (i + 1 < TAM_POPULACAO)
                    nova_geracao[i + 1] = filho2;
            }
            for (int i = 0; i < TAM_POPULACAO; i++)
                populacao[i] = nova_geracao[i];
            aplicar_morte(populacao, geracao);
            // Melhor local
            melhor_local = populacao[0];
            for (int i = 1; i < TAM_POPULACAO; i++)
                if (populacao[i].fitness > melhor_local.fitness)
                    melhor_local = populacao[i];
            imprimir_populacao(populacao, geracao, rank-1);
            // Migração
            if (geracao % INTERVALO_MIGRACAO == 0) {
                // Envia melhor para master
                MPI_Send(&melhor_local, 1, MPI_INDIVIDUO, 0, 0, MPI_COMM_WORLD);
                // Recebe migrante do master
                MPI_Recv(&migrante, 1, MPI_INDIVIDUO, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // Substitui pior
                int pior = 0;
                for (int i = 1; i < TAM_POPULACAO; i++)
                    if (populacao[i].fitness < populacao[pior].fitness)
                        pior = i;
                printf("[ILHA %d] Substituindo indivíduo com fitness %.2f por migrante com fitness %.2f\n",
                    rank, populacao[pior].fitness, migrante.fitness);
                populacao[pior] = migrante;
            }
            // Log
            if (geracao % 5 == 0) {
                printf("\n==== Relatório Geração %d (Ilha %d) ====
Melhor Fitness: %.2f\n", geracao, rank, melhor_local.fitness);
            }
        }
    }
    MPI_Type_free(&MPI_INDIVIDUO);
    MPI_Finalize();
    return 0;
}

// Funções auxiliares (iguais ao cod_seq.c)
int binario_para_inteiro(int genes[], int inicio, int fim) {
    int valor = 0;
    for (int i = inicio; i <= fim; i++) {
        valor = (valor << 1) | genes[i];
    }
    return valor;
}

float calcular_fitness(Individuo individuo) {
    int rad = binario_para_inteiro(individuo.genes, 0, 7);
    int fome = binario_para_inteiro(individuo.genes, 8, 15);
    int energia = binario_para_inteiro(individuo.genes, 16, 23);
    int sanidade = binario_para_inteiro(individuo.genes, 24, 31);
    int doenca = binario_para_inteiro(individuo.genes, 32, 39);
    float somaPesos = PESO_RAD + PESO_FOME + PESO_ENERGIA + PESO_SANIDADE + PESO_DOENCA;
    float fitness = (
        ((255 - rad) / 255.0) * PESO_RAD +
        ((255 - fome) / 255.0) * PESO_FOME +
        (energia / 255.0) * PESO_ENERGIA +
        (sanidade / 255.0) * PESO_SANIDADE +
        ((255 - doenca) / 255.0) * PESO_DOENCA
    ) * (100 / somaPesos);
    return fitness;
}

void inicializar_populacao(Individuo populacao[]) {
    for (int i = 0; i < TAM_POPULACAO; i++) {
        for (int j = 0; j < TAM_CROMOSSOMO; j++) {
            populacao[i].genes[j] = rand() % 2;
        }
        populacao[i].fitness = calcular_fitness(populacao[i]);
    }
}

Individuo torneio(Individuo populacao[]) {
    int a = rand() % TAM_POPULACAO;
    int b = rand() % TAM_POPULACAO;
    return (populacao[a].fitness > populacao[b].fitness) ? populacao[a] : populacao[b];
}

void cruzamento(Individuo pai, Individuo mae, Individuo *filho1, Individuo *filho2) {
    if ((float)rand() / RAND_MAX < PROB_CROSSOVER) {
        int ponto = rand() % TAM_CROMOSSOMO;
        for (int i = 0; i < ponto; i++) {
            filho1->genes[i] = pai.genes[i];
            filho2->genes[i] = mae.genes[i];
        }
        for (int i = ponto; i < TAM_CROMOSSOMO; i++) {
            filho1->genes[i] = mae.genes[i];
            filho2->genes[i] = pai.genes[i];
        }
    } else {
        *filho1 = pai;
        *filho2 = mae;
    }
    filho1->fitness = calcular_fitness(*filho1);
    filho2->fitness = calcular_fitness(*filho2);
}

Individuo mutacao(Individuo individuo) {
    for (int i = 0; i < TAM_CROMOSSOMO; i++) {
        if ((float)rand() / RAND_MAX < PROB_MUTACAO) {
            individuo.genes[i] = !individuo.genes[i];
        }
    }
    individuo.fitness = calcular_fitness(individuo);
    return individuo;
}

void imprimir_individuo(Individuo ind) {
    printf("Genes: ");
    for (int i = 0; i < TAM_CROMOSSOMO; i++) {
        printf("%d", ind.genes[i]);
    }
    printf("\n");
    int rad = binario_para_inteiro(ind.genes, 0, 7);
    int fome = binario_para_inteiro(ind.genes, 8, 15);
    int energia = binario_para_inteiro(ind.genes, 16, 23);
    int sanidade = binario_para_inteiro(ind.genes, 24, 31);
    int doenca = binario_para_inteiro(ind.genes, 32, 39);
    int rad_pct = rad / 255.0 * 100;
    int fome_pct = fome / 255.0 * 100;
    int energia_pct = energia / 255.0 * 100;
    int sanidade_pct = sanidade / 255.0 * 100;
    int doenca_pct = doenca / 255.0 * 100;
    printf("Radiacao: %d%%, Fome: %d%%, Energia: %d%%, Sanidade: %d%%, Doenca: %d%% | Fitness: %.2f\n",
       rad_pct, fome_pct, energia_pct, sanidade_pct, doenca_pct, ind.fitness);
}

void imprimir_populacao(Individuo populacao[], int geracao, int ilha) {
    printf("\n==== População da Ilha %d - Geração %d ====\n", ilha + 1, geracao);
    for (int i = 0; i < TAM_POPULACAO; i++) {
        printf("Indivíduo %d: ", i);
        imprimir_individuo(populacao[i]);
    }
}

void aplicar_morte(Individuo populacao[], int geracao) {
    for (int i = 0; i < TAM_POPULACAO; i++) {
        if (populacao[i].fitness < LIMIAR_MORTE) {
            printf("\n[Geracao %d] Individuo %d morreu (fitness %.2f abaixo do limiar %d)\n", 
                    geracao, i, populacao[i].fitness, LIMIAR_MORTE);
            printf("[Geracao %d] Nova expedicao buscando novo individuo %d com genes: ", geracao, i);
            for (int j = 0; j < TAM_CROMOSSOMO; j++) {
                populacao[i].genes[j] = rand() % 2;
                printf("%d", populacao[i].genes[j]);
            }
            populacao[i].fitness = calcular_fitness(populacao[i]);
            printf("\n[Geracao %d] Novo individuo %d resgatado pela colonia com fitness %.2f\n", 
                    geracao, i, populacao[i].fitness);
        }
    }
}

void inicializar_parametros() {
    PROB_CROSSOVER = (float)rand() / RAND_MAX;
    PROB_MUTACAO = (float)rand() / RAND_MAX;
    LIMIAR_MORTE = rand() % 101;
    float pesos[5], soma = 0.0;
    for (int i = 0; i < 5; i++) {
        pesos[i] = (float)rand() / RAND_MAX;
        soma += pesos[i];
    }
    PESO_RAD = pesos[0] / soma;
    PESO_FOME = pesos[1] / soma;
    PESO_ENERGIA = pesos[2] / soma;
    PESO_SANIDADE = pesos[3] / soma;
    PESO_DOENCA = pesos[4] / soma;
    printf("\n==== Estatistica - Dadas ====\n");
    printf("CHEGADAS DE DADOS PARA NOVA POPULACAO\n");
    printf("PROB_CROSSOVER: %.2f\n", PROB_CROSSOVER);
    printf("PROB_MUTACAO: %.2f\n", PROB_MUTACAO);
    printf("LIMIAR_MORTE: %d\n", LIMIAR_MORTE);
    printf("PESOS: RAD=%.2f, FOME=%.2f, ENERGIA=%.2f, SANIDADE=%.2f, DOENCA=%.2f (soma=%.2f)\n",
        PESO_RAD, PESO_FOME, PESO_ENERGIA, PESO_SANIDADE, PESO_DOENCA,
        PESO_RAD + PESO_FOME + PESO_ENERGIA + PESO_SANIDADE + PESO_DOENCA);
} 