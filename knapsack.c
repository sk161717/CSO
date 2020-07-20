#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

// 個数の上限はあらかじめ定めておく
const int max_items = 1000;
const double MR = 0.2;  //mixture ratio of tracing cat
const int SMP = 3;      //SMP indicate how many points is explored
const double CDC = 0.2; //indicate how many dimension is varied in seeking mode
const int SPC = 1;      //decide whether current position can be one of candidates
const double PMO = 0.2; //indicate the probability of mutation on the selected dimension
const double w = 0.5;
const double c1 = 1;
const double vmax = 4;
const int rand_seed = 9;

// 以下は構造体の定義と関数のプロトタイプ宣言

// 構造体 itemset
// number分の価値value と 重さweight を格納しておく
// データをポインタで定義しており、mallocする必要あり
typedef struct itemset
{
  int number;
  double *value;
  double *weight;
} Itemset;

typedef struct cat
{
  int *flags;
  double *vkd1;
  double *vkd0;
  double *vkd;
  int state; //if state==1 trace else seek
  double fitness;
} Cat;

typedef struct cats
{
  Cat *cats;
} Cats;

// 関数のプロトサイプ宣言

// Itemset *init_itemset(int, int);
//
// itemsetを初期化し、そのポインタを返す関数
// 引数:
//  品物の個数: number (int)
//  乱数シード: seed (int) // 品物の値をランダムにする
// 返り値:
//  確保されたItemset へのポインタ
Itemset *init_itemset(int number);

// void free_itemset();

// Itemset *load_itemset(char *filename)
//
// ファイルからItemset を設定し、確保された領域へのポインタを返す関数 [未実装]
// 引数:
//  Itemsetの必要パラメータが記述されたバイナリファイルのファイル名 filename (char*)
// 返り値:
//  Itemset へのポインタ
Itemset *load_itemset(char *filename);

// void print_itemset(const Itemset *list)
//
// Itemsetの内容を標準出力に表示する関数
void print_itemset(const Itemset *list);

// void save_itemset(char *filename)
//
// Itemsetのパラメータを記録したバイナリファイルを出力する関数
// 引数:
// Itemsetの必要パラメータを吐き出すファイルの名前 filename (char*)
// 返り値:
//  なし
void save_itemset(char *filename);

// main関数
// プログラム使用例: ./knapsack 10 20
//  10個の品物を設定し、キャパ20 でナップサック問題をとく
Cats *initialize_cats(int cat_size, int dimension);

double rand_01()
{
  return rand() / (RAND_MAX * 1.0);
}

double sigmoid(double x)
{
  return 1.0 / (1.0 + exp((-1.0) * x));
}

double calculate_fitness(Cat *cat, int dimension, Itemset *items, double W)
{
  double total_value = 0;
  double total_weight = 0;
  for (int i = 0; i < dimension; i++)
  {
    total_value += cat->flags[i] * items->value[i];
    total_weight += cat->flags[i] * items->weight[i];
  }
  if (total_weight > W)
  {
    total_value = 0;
  }
  cat->fitness = total_value;
  return total_value;
}

Cat *bestfitness_cat(Cats *cats, int dimension, Itemset *items, double W, int cat_size)
{
  Cat *best_cat;
  Cat *best_cat_cp;
  best_cat_cp = (Cat *)malloc(sizeof(Cat));
  best_cat_cp->flags = (int *)malloc(sizeof(int) * dimension);

  double best_cat_fitness = -1;
  for (int i = 0; i < cat_size; i++)
  {
    double tmp_fitness = calculate_fitness(&(cats->cats[i]), dimension, items, W);
    if (tmp_fitness > best_cat_fitness)
    {
      best_cat = &(cats->cats[i]);
      best_cat_fitness = tmp_fitness;
    }
  }
  for (int i = 0; i < dimension; i++)
  {
    best_cat_cp->flags[i] = best_cat->flags[i];
  }
  return best_cat_cp;
}

double calculate_t(Cat *best_fit_cat, Cat *cat, int d)
{
  double dkd1, dkd0;
  double r1c1 = c1 * rand() / (RAND_MAX * 1.0);
  if (best_fit_cat->flags[d] == 1)
  {
    dkd1 = r1c1;
    dkd0 = -r1c1;
  }
  else
  {
    dkd1 = -r1c1;
    dkd0 = r1c1;
  }
  cat->vkd1[d] = w * cat->vkd1[d] + dkd1;
  cat->vkd0[d] = w * cat->vkd0[d] + dkd0;
  if (fabs(cat->vkd[d]) > vmax)
  {
    cat->vkd1[d] = vmax * cat->vkd1[d] / fabs(cat->vkd1[d]);
    cat->vkd0[d] = vmax * cat->vkd0[d] / fabs(cat->vkd0[d]);
  }
  if (cat->flags[d] == 0)
  {
    cat->vkd[d] = cat->vkd1[d];
  }
  else
  {
    cat->vkd[d] = cat->vkd0[d];
  }
  return sigmoid(cat->vkd[d]);
}

int *randsample(int dimension, int variable_dimension)
{
  int *rand_array = (int *)calloc(dimension, sizeof(int));

  while (1)
  {
    int random = (int)round(rand() / (RAND_MAX * 1.0) * dimension);
    if (rand_array[random] == 0)
    {
      rand_array[random] = 1;
      variable_dimension--;
    }
    if (variable_dimension == 0)
    {
      break;
    }
  }
  return rand_array;
}

int roulette(double *probability)
{
  double sum;
  for (int i = 0; i < SMP; i++)
  {
    sum += probability[i];
  }
  double random = rand_01() * sum;
  double tmp = 0;
  for (int i = 0; i < SMP; i++)
  {
    tmp += probability[i];
    if (random < tmp)
    {
      return i;
    }
  }
}

void seek(Cat *cat, int dimension, Itemset *items, double W)
{
  Cats *copy_cats = (Cats *)malloc(sizeof(Cats));
  copy_cats->cats = (Cat *)malloc(sizeof(Cat) * SMP);
  for (int i = 0; i < SMP; i++)
  {
    copy_cats->cats[i].flags = (int *)malloc(sizeof(int) * dimension);
    for (int j = 0; j < dimension; j++)
    {
      copy_cats->cats[i].flags[j] = cat->flags[j];
    }
  }
  int variable_dimension = (int)round(CDC * dimension);
  int variable_cats = SMP - SPC;

  for (int i = 0; i < variable_cats; i++)
  {
    int *rand_array = randsample(dimension, variable_dimension);
    for (int j = 0; j < dimension; j++)
    {
      if (rand_array[j] == 1)
      {
        double rand = rand_01();
        if (rand < PMO)
        {
          copy_cats->cats[i].flags[j] = 1 - copy_cats->cats[i].flags[j];
        }
      }
    }
  }
  double FSMax = 0, FSMin = RAND_MAX;
  for (int i = 0; i < SMP; i++)
  {
    copy_cats->cats[i].fitness = calculate_fitness(&(copy_cats->cats[i]), dimension, items, W);
    if (copy_cats->cats[i].fitness > FSMax)
    {
      FSMax = copy_cats->cats[i].fitness;
    }
    if (copy_cats->cats[i].fitness < FSMin)
    {
      FSMin = copy_cats->cats[i].fitness;
    }
  }
  double *probability = (double *)malloc(sizeof(double) * SMP);
  if (FSMax == FSMin)
  {
    for (size_t i = 0; i < SMP; i++)
    {
      probability[i] = 1;
    }
  }
  else
  {
    for (size_t i = 0; i < SMP; i++)
    {
      probability[i] = (copy_cats->cats[i].fitness - FSMin) / (FSMax - FSMin);
    }
  }
  int index = roulette(probability);
  for (int i = 0; i < dimension; i++)
  {
    cat->flags[i] = copy_cats->cats[index].flags[i];
  }
}

void move_phase(Cats *cats, Cat *best_fit_cat, int cat_size, int dimension, Itemset *items, double W)
{
  for (size_t i = 0; i < cat_size; i++)
  {
    if (cats->cats[i].state == 1)
    {
      for (size_t j = 0; j < dimension; j++)
      {

        double random = rand() / (RAND_MAX * 1.0);
        double t = calculate_t(best_fit_cat, &(cats->cats[i]), j);
        if (t > random)
        {
          cats->cats[i].flags[j] = best_fit_cat->flags[j];
        }
        //else no change
      }
    }
    else
    {
      seek(&(cats->cats[i]), dimension, items, W);
    }
  }
}

void print_status(Cats *cats, int dimension, int cat_size)
{
  for (int i = 0; i < cat_size; i++)
  {
    for (int j = 0; j < dimension; j++)
    {
      printf("%d", cats->cats[i].flags[j]);
    }
    /*for (int j = 0; j < dimension; j++)
    {
      printf(" %lf", cats->cats[i].vkd0[j]);
    }
    for (int j = 0; j < dimension; j++)
    {
      printf(" %lf", cats->cats[i].vkd1[j]);
    }*/
    printf("state is %d", cats->cats[i].state);
    printf(" %lf", cats->cats[i].fitness);
    printf("\n");
  }
}

void state_change(Cats *cats, int dimension, int cat_size)
{
  for (int i = 0; i < cat_size; i++)
  {
    double random_MR = rand() / (RAND_MAX * 1.0);
    if (random_MR < MR)
    {
      cats->cats[i].state = 1;
    }
    else
    {
      cats->cats[i].state = 0;
    }
  }
}

int main(int argc, char **argv)
{
  clock_t start, end;

  /* 引数処理: ユーザ入力が正しくない場合は使い方を標準エラーに表示して終了 */
  if (argc != 3)
  {
    fprintf(stderr, "usage: %s <the number of items (int)> <max capacity (double)>\n", argv[0]);
    exit(1);
  }

  const int dimension = atoi(argv[1]);
  assert(dimension <= max_items);

  const double W = atof(argv[2]);
  int rand_seed[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  printf("max capacity: W = %.f, # of items: %d\n", W, dimension);
  Itemset *items = init_itemset(dimension);
  int cat_size = dimension * 40;
  for (int a = 0; a < 10; a++)
  {
    srand(rand_seed[a]);
    start = clock();

    //print_itemset(items);

    Cats *cats = initialize_cats(cat_size, dimension);
    Cat *best_cat = bestfitness_cat(cats, dimension, items, W, cat_size);
    printf("0= %lf\n", calculate_fitness(best_cat, dimension, items, W));
    //print_status(cats, dimension, cat_size);
    while (1)
    {

      best_cat = bestfitness_cat(cats, dimension, items, W, cat_size);
      move_phase(cats, best_cat, cat_size, dimension, items, W);
      state_change(cats, dimension, cat_size);
      //printf("%lf\n", calculate_fitness(best_cat, dimension, items, W));

      best_cat = bestfitness_cat(cats, dimension, items, W, cat_size);
      if (calculate_fitness(best_cat, dimension, items, W) == 13549094)
      {
        break;
      }
    }

    best_cat = bestfitness_cat(cats, dimension, items, W, cat_size);
    printf("\n");
    //print_status(cats, dimension, cat_size);
    printf("%lf\n", calculate_fitness(best_cat, dimension, items, W));
    end = clock();
    printf("%.2f秒かかりました\n", (double)(end - start) / CLOCKS_PER_SEC);
    free(cats);
    free(best_cat);
  }

  return 0;
}

Cats *initialize_cats(int cat_size, int dimension)
{
  Cats *cats;
  cats = (Cats *)malloc(sizeof(Cats));
  cats->cats = (Cat *)malloc(sizeof(Cat) * cat_size);
  for (int i = 0; i < cat_size; i++)
  {
    cats->cats[i].flags = (int *)malloc(sizeof(int) * dimension);
    cats->cats[i].vkd1 = (double *)malloc(sizeof(double) * dimension);
    cats->cats[i].vkd0 = (double *)malloc(sizeof(double) * dimension);
    cats->cats[i].vkd = (double *)malloc(sizeof(double) * dimension);
    for (int j = 0; j < dimension; j++)
    {
      double random_vkd1 = rand() / (RAND_MAX * 1.0);
      cats->cats[i].vkd1[j] = random_vkd1 - 0.5;
      double random_vkd0 = rand() / (RAND_MAX * 1.0);
      cats->cats[i].vkd0[j] = random_vkd0 - 0.5;
      double random_vkd = rand() / (RAND_MAX * 1.0);
      cats->cats[i].vkd[j] = random_vkd - 0.5;

      double random = rand() / (RAND_MAX * 1.0);
      int flag;
      if (random > 0.5)
      {
        flag = 1;
      }
      else
      {
        flag = 0;
      }

      cats->cats[i].flags[j] = flag;
    }

    double random_MR = rand() / (RAND_MAX * 1.0);
    if (random_MR < MR)
    {
      cats->cats[i].state = 1;
    }
    else
    {
      cats->cats[i].state = 0;
    }
  }
  return cats;
}

// 構造体をポインタで確保するお作法を確認してみよう
Itemset *init_itemset(int number)
{
  Itemset *list = (Itemset *)malloc(sizeof(Itemset));

  list->number = number;

  list->value = (double *)malloc(sizeof(double) * number);
  list->weight = (double *)malloc(sizeof(double) * number);

  FILE *fp_value;
  FILE *fp_weight;

  if ((fp_value = fopen("p08_p.txt", "r")) == NULL)
  {
    printf("p07_p.txt error");
    exit(1);
  }
  if ((fp_weight = fopen("p08_w.txt", "r")) == NULL)
  {
    printf("p07_w.txt error");
    exit(1);
  };

  for (int i = 0; i < number; i++)
  {
    double *tmp_v = (double *)malloc(sizeof(double));
    double *tmp_w = (double *)malloc(sizeof(double));
    fscanf(fp_value, "%lf", tmp_v);
    fscanf(fp_weight, "%lf", tmp_w);
    list->value[i] = *tmp_v;
    list->weight[i] = *tmp_w;
  }
  fclose(fp_value);
  fclose(fp_weight);

  return list;
}

// itemset の free関数
void free_itemset(Itemset *list)
{
  free(list->value);
  free(list->weight);
  free(list);
}

// 表示関数
void print_itemset(const Itemset *list)
{
  int n = list->number;
  for (int i = 0; i < n; i++)
  {
    printf("v[%d] = %4.1f, w[%d] = %4.1f\n", i, list->value[i], i, list->weight[i]);
  }
  printf("----\n");
}
