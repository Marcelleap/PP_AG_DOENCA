#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define BITS_POR_GENE 8
#define NUM_GENES 5
#define TAM_POPULACAO 5
#define TAM_CROMOSSOMO (BITS_POR_GENE * NUM_GENES)
#define MAX_GER 10
#define RAD_MAX 2

#define NUM_ILHAS 3
#define INTERVALO_MIGRACAO 3

float PROB_CROSSOVER;
float PROB_MUTACAO;
int LIMIAR_MORTE;
float PESO_RAD, PESO_FOME, PESO_ENERGIA, PESO_SANIDADE, PESO_DOENCA;

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

int main() {
    srand(time(NULL));
    inicializar_parametros();

    clock_t start_time = clock(); // INÍCIO DA CONTAGEM

    Individuo ilhas[NUM_ILHAS][TAM_POPULACAO];
    Individuo nova_geracao[NUM_ILHAS][TAM_POPULACAO];
    Individuo melhores_ilhas[NUM_ILHAS];
    Individuo migrantes[NUM_ILHAS];

    // Inicializar todas as ilhas
    for (int ilha = 0; ilha < NUM_ILHAS; ilha++)
        inicializar_populacao(ilhas[ilha]);

    for (int geracao = 0; geracao < MAX_GER; geracao++) {
        // Evolução local de cada ilha
        for (int ilha = 0; ilha < NUM_ILHAS; ilha++) {
            for (int i = 0; i < TAM_POPULACAO; i += 2) {
                Individuo pai = torneio(ilhas[ilha]);
                Individuo mae = torneio(ilhas[ilha]);
                Individuo filho1, filho2;
                cruzamento(pai, mae, &filho1, &filho2);
                filho1 = mutacao(filho1);
                filho2 = mutacao(filho2);
                nova_geracao[ilha][i] = filho1;
                if (i + 1 < TAM_POPULACAO)
                    nova_geracao[ilha][i + 1] = filho2;
            }
            for (int i = 0; i < TAM_POPULACAO; i++)
                ilhas[ilha][i] = nova_geracao[ilha][i];
            aplicar_morte(ilhas[ilha], geracao);
            // Encontrar melhor local
            melhores_ilhas[ilha] = ilhas[ilha][0];
            for (int i = 1; i < TAM_POPULACAO; i++)
                if (ilhas[ilha][i].fitness > melhores_ilhas[ilha].fitness)
                    melhores_ilhas[ilha] = ilhas[ilha][i];

            // Imprimir população da ilha após evolução
            imprimir_populacao(ilhas[ilha], geracao, ilha);
        }

         // --- Nó zero imprime os melhores de cada ilha ---
        printf("\n[NO ZERO] Melhores indivíduos recebidos na geração %d:\n", geracao);
        for (int ilha = 0; ilha < NUM_ILHAS; ilha++) {
            printf("Ilha %d: ", ilha + 1);
            imprimir_individuo(melhores_ilhas[ilha]);
        }

        // --- Migração circular a cada INTERVALO_MIGRACAO ---
        if (geracao % INTERVALO_MIGRACAO == 0) {
            // Ordenar índices das ilhas por fitness (decrescente)
            int indices[NUM_ILHAS] = {0, 1, 2};
            for (int i = 0; i < NUM_ILHAS - 1; i++) {
                for (int j = i + 1; j < NUM_ILHAS; j++) {
                    if (melhores_ilhas[indices[j]].fitness > melhores_ilhas[indices[i]].fitness) {
                        int tmp = indices[i];
                        indices[i] = indices[j];
                        indices[j] = tmp;
                    }
                }
            }
            migrantes[0] = melhores_ilhas[indices[2]]; // 3º melhor vai para ilha 0
            migrantes[1] = melhores_ilhas[indices[0]]; // 1º melhor vai para ilha 1
            migrantes[2] = melhores_ilhas[indices[1]]; // 2º melhor vai para ilha 2

            // Substituir o pior indivíduo de cada ilha pelo migrante
            for (int ilha = 0; ilha < NUM_ILHAS; ilha++) {
                int pior = 0;
                for (int i = 1; i < TAM_POPULACAO; i++)
                    if (ilhas[ilha][i].fitness < ilhas[ilha][pior].fitness)
                        pior = i;
                printf("[ILHA %d] Substituindo indivíduo com fitness %.2f por migrante com fitness %.2f\n",
                    ilha + 1, ilhas[ilha][pior].fitness, migrantes[ilha].fitness);
                ilhas[ilha][pior] = migrantes[ilha];
            }
        }

        // Log de cada ilha
        if (geracao % 5 == 0) {
            printf("\n==== Relatório Geração %d ====\n", geracao);
            for (int ilha = 0; ilha < NUM_ILHAS; ilha++) {
                printf("Ilha %d - Melhor Fitness: %.2f\n", ilha + 1, melhores_ilhas[ilha].fitness);
            }
        }
    }

    clock_t end_time = clock();
    double tempo_execucao = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("\nTEMPO DE EXECUÇÃO SEQUENCIAL: %.6f segundos\n", tempo_execucao);

    return 0;
}

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
