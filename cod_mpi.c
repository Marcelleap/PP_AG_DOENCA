#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"


#define BITS_POR_GENE 8
#define NUM_GENES 5
#define TAM_POPULACAO 5
#define TAM_CROMOSSOMO (BITS_POR_GENE * NUM_GENES)
#define MAX_GER 10
#define RAD_MAX 2


float PROB_CROSSOVER ;
float PROB_MUTACAO ;

int LIMIAR_MORTE ; 

float PESO_RAD;
float PESO_FOME;
float PESO_ENERGIA;
float PESO_SANIDADE;
float PESO_DOENCA;
    

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
void imprimir_populacao (Individuo populacao[], int geracao);
void imprimir_individuo(Individuo ind); 
void aplicar_morte(Individuo populacao[], int geracao);
void inicializar_parametros();
int ilha_processo(int rank);
void master_controlador();


// Definições existentes...
#define TAG_MELHORES 1
#define TAG_MIGRACAO 2
#define INTERVALO_MIGRACAO 3  // Número de gerações entre migrações

int main(int argc, char *argv[]) {
    int rank, size;
	double start_time, end_time, resolution;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 4) {
        if (rank == 0) {
            fprintf(stderr, "Este programa deve ser executado com exatamente 4 processos!\n");
        }
        MPI_Finalize();
        return 1;
    }

	// Obter a resolucao do relogio
    resolution = MPI_Wtick();
    if (rank == 0) {
        printf("Resolucao do relogio MPI: %.9f segundos\n", resolution);
        start_time = MPI_Wtime();  // Inicia cronometro apenas no master
    }

    srand(time(NULL) + rank);  // Semente diferente para cada nó
    inicializar_parametros();

    if (rank == 0) {
		master_controlador();
    } else {
        ilha_processo(rank);
    }

	if (rank == 0) {
        end_time = MPI_Wtime();
        double tempo_paralelo = end_time - start_time;
        printf("\nTEMPO DE EXECUÇÃO PARALELO: %.6f segundos (resolução: %.9f s)\n", 
              tempo_paralelo, resolution);
        
        FILE *f = fopen("tempos.txt", "a");
        fprintf(f, "PARALELO,%.6f,%.9f\n", tempo_paralelo, resolution);
        fclose(f);
    }


    MPI_Finalize();
    return 0;
}

void master_controlador() {
    printf("=== MASTER INICIANDO CONTROLE DE MIGRAÇÃO ENTRE ILHAS ===\n");
    
    Individuo melhores_ilhas[3];  // Armazena os melhores de cada ilha
    MPI_Request requests[3];
    int ilhas_ativas = 3;
    
    for (int geracao = 0; geracao < MAX_GER; geracao++) {
        // 1. Receber os melhores indivíduos de cada ilha (não bloqueante)
        for (int i = 0; i < 3; i++) {
            MPI_Irecv(&melhores_ilhas[i], sizeof(Individuo), MPI_BYTE, 
                     i+1, TAG_MELHORES, MPI_COMM_WORLD, &requests[i]);
        }
        
            // Imprimir população da ilha após evolução
        //imprimir_populacao(ilhas[ilha], geracao, ilha);

        // 2. Processar enquanto aguarda (pode fazer análise ou logging)
        MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);

        printf("\n[MASTER] Melhores indivíduos recebidos na geração %d:\n", geracao);
        for (int i = 0; i < 3; i++) {
            printf("Ilha %d: ", i+1);
            imprimir_individuo(melhores_ilhas[i]);
        }
        
        // 3. Selecionar os melhores para migração (estratégia elitista)
        Individuo melhores_para_migrar[3];
        int indices[3] = {0, 1, 2};
        
        // Ordenar ilhas por fitness do melhor indivíduo
        for (int i = 0; i < 2; i++) {
            for (int j = i+1; j < 3; j++) {
                if (melhores_ilhas[indices[j]].fitness > melhores_ilhas[indices[i]].fitness) {
                    int temp = indices[i];
                    indices[i] = indices[j];
                    indices[j] = temp;
                }
            }
        }
        
        // 4. Preparar migração circular: 1->2->3->1
        if (geracao % INTERVALO_MIGRACAO == 0) {
            printf("\n[MASTER] Coordenando migração na geração %d\n", geracao);
            
            melhores_para_migrar[0] = melhores_ilhas[indices[0]];  // Melhor absoluto
            melhores_para_migrar[1] = melhores_ilhas[indices[1]];  // Segundo melhor
            melhores_para_migrar[2] = melhores_ilhas[indices[2]];  // Terceiro melhor
            
            // Enviar indivíduos para migração controlada
            MPI_Send(&melhores_para_migrar[2], sizeof(Individuo), MPI_BYTE, 1, TAG_MIGRACAO, MPI_COMM_WORLD);  // 3º melhor vai para ilha 1
            MPI_Send(&melhores_para_migrar[0], sizeof(Individuo), MPI_BYTE, 2, TAG_MIGRACAO, MPI_COMM_WORLD);  // 1º melhor vai para ilha 2
            MPI_Send(&melhores_para_migrar[1], sizeof(Individuo), MPI_BYTE, 3, TAG_MIGRACAO, MPI_COMM_WORLD);  // 2º melhor vai para ilha 3
            
            printf("[MASTER] Migração completada:\n");
            printf(" - Ilha 1 recebeu indivíduo com fitness %.2f\n", melhores_para_migrar[2].fitness);
            printf(" - Ilha 2 recebeu indivíduo com fitness %.2f\n", melhores_para_migrar[0].fitness);
            printf(" - Ilha 3 recebeu indivíduo com fitness %.2f\n", melhores_para_migrar[1].fitness);
        }
        
        // 5. Exibir relatório periódico
        if (geracao % 5 == 0) {
            printf("\n[MASTER] Relatório Geração %d:\n", geracao);
            for (int i = 0; i < 3; i++) {
                printf("Ilha %d - Melhor Fitness: %.2f\n", i+1, melhores_ilhas[i].fitness);
            }
        }
    }
}

int ilha_processo(int rank) {
    Individuo populacao[TAM_POPULACAO];
    Individuo nova_geracao[TAM_POPULACAO];
    Individuo melhor_local;
    Individuo migrante;
    
    inicializar_populacao(populacao);
    
    for (int geracao = 0; geracao < MAX_GER; geracao++) {
        // 1. Evolução local
        for (int i = 0; i < TAM_POPULACAO; i += 2) {
            Individuo pai = torneio(populacao);
            Individuo mae = torneio(populacao);
            
            Individuo filho1, filho2;
            cruzamento(pai, mae, &filho1, &filho2);
            
            filho1 = mutacao(filho1);
            filho2 = mutacao(filho2);
            
            nova_geracao[i] = filho1;
            nova_geracao[i + 1] = filho2;
        }
        
        // 2. Substituição geracional
        for (int i = 0; i < TAM_POPULACAO; i++) {
            populacao[i] = nova_geracao[i];
        }
        
        aplicar_morte(populacao, geracao);
        
        // 3. Encontrar o melhor local para enviar ao master
        melhor_local = populacao[0];
        for (int i = 1; i < TAM_POPULACAO; i++) {
            if (populacao[i].fitness > melhor_local.fitness) {
                melhor_local = populacao[i];
            }
        }
        MPI_Send(&melhor_local, sizeof(Individuo), MPI_BYTE, 0, TAG_MELHORES, MPI_COMM_WORLD);
        
        // 4. Receber migrantes do master (se for época de migração)
        if (geracao % INTERVALO_MIGRACAO == 0) {
            MPI_Recv(&migrante, sizeof(Individuo), MPI_BYTE, 0, TAG_MIGRACAO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Substituir o pior indivíduo
            int pior = 0;
            for (int i = 1; i < TAM_POPULACAO; i++) {
                if (populacao[i].fitness < populacao[pior].fitness) {
                    pior = i;
                }
            }
            
            printf("[ILHA %d] Substituindo indivíduo com fitness %.2f por migrante com fitness %.2f\n",
                  rank, populacao[pior].fitness, migrante.fitness);
            
            populacao[pior] = migrante;
        }
        
        // 5. Log local (opcional)
        if (geracao % 10 == 0) {
            printf("[ILHA %d] Geração %d - Melhor local: %.2f\n",
                  rank, geracao, melhor_local.fitness);
        }
    }
}


// Função para converter binário em inteiro
int binario_para_inteiro(int genes[], int inicio, int fim) {
    
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

// Inicializar população
void inicializar_populacao(Individuo indviduo[]) {
    
    for (int i = 0; i < TAM_POPULACAO; i++) {
        
        for (int j = 0; j < TAM_CROMOSSOMO; j++) {
            indviduo[i].genes[j] = rand() % 2;
        }
        indviduo[i].fitness = calcular_fitness(indviduo[i]);
    }

}

void imprimir_populacao(Individuo populacao[], int geracao, int ilha) {
    printf("\n==== População da Ilha %d - Geração %d ====\n", ilha + 1, geracao);
    for (int i = 0; i < TAM_POPULACAO; i++) {
        printf("Indivíduo %d: ", i);
        imprimir_individuo(populacao[i]);
    }
}

// Seleção por torneio
Individuo torneio(Individuo populacao[]) {
    
    int a = rand() % TAM_POPULACAO;
    int b = rand() % TAM_POPULACAO;

    if (populacao[a].fitness > populacao[b].fitness) {
        return populacao[a];
    } else {
        return populacao[b];
    }
}

// Cruzamento de um ponto
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
    } 
    
    else {
        *filho1 = pai;
        *filho2 = mae;
    }

    filho1->fitness = calcular_fitness(*filho1);
    filho2->fitness = calcular_fitness(*filho2);
}

// Mutação
Individuo mutacao(Individuo individuo) {
    
    for (int i = 0; i < TAM_CROMOSSOMO; i++) {
        if ((float)rand() / RAND_MAX < PROB_MUTACAO) {
            individuo.genes[i] = !individuo.genes[i];
        }
    }
    individuo.fitness = calcular_fitness(individuo);
    
    return individuo;
}
 /*
void imprimir_populacao (Individuo populacao[], int geracao)
{
        printf("\n==== Populacao - Geracao %d ====\n", geracao);
        
        for (int i = 0; i<TAM_POPULACAO; i++){
            printf("Individuo %d\n", i);
            imprimir_individuo(populacao[i]); 
            printf("\n");
        }
}
*/

// Impressão/
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
    // Probabilidades entre 0 e 1
    PROB_CROSSOVER = (float)rand() / RAND_MAX;
    PROB_MUTACAO = (float)rand() / RAND_MAX;

    // Limiar de morte entre 0 e 100
    LIMIAR_MORTE = rand() % 101;

    // Pesos aleatórios e normalização
    float pesos[5];
    float soma = 0.0;
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
    // Exibir os valores sorteados para conferência
    printf("CHEGADAS DE DADOS PARA NOVA POPULACAO\n");
    printf("PROB_CROSSOVER: %.2f\n", PROB_CROSSOVER);
    printf("PROB_MUTACAO: %.2f\n", PROB_MUTACAO);
    printf("LIMIAR_MORTE: %d\n", LIMIAR_MORTE);
    printf("PESOS: RAD=%.2f, FOME=%.2f, ENERGIA=%.2f, SANIDADE=%.2f, DOENCA=%.2f (soma=%.2f)\n",
        PESO_RAD, PESO_FOME, PESO_ENERGIA, PESO_SANIDADE, PESO_DOENCA,
        PESO_RAD + PESO_FOME + PESO_ENERGIA + PESO_SANIDADE + PESO_DOENCA);
}
