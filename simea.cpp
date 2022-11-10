#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "header4.h"
#include "paramt3.h"




/*軌道計算用の関数の練習*/

/*double GG_K(double k,double et, double t_ini){
        double fk;*/
/*if((et-t_ini)<=1827){
    fk=0;
}else{*/

    /*fk=1/(1+k);
//}
return fk;
}*/
double GG_K(double k){
        double fk;
/*if((et-t_ini)<=1096){
    fk=0;
}else{*/

    fk=1/(1+fabs(k));
//}
return fk;
}

/*メモリ確保のための関数。確保できなかったら場合はエラー・終了*/
void* my_malloc(int size){
    void* ptr = malloc(size);
    if(!ptr){
        fprintf(stderr,"failed to malloc.\n");
        exit (1);
    }
    return ptr;
}

/*gtype(遺伝子コード)の実数GAのための実装*/
/*gtype(int配列)を作り最初のアドレスを返す*/
gtype_t mk_gtype(int code_length){
    gtype_t gtype=(gtype_t)my_malloc(sizeof(int)*code_length);
    return gtype;
}

/*gtypeのメモリを開放する*/
void free_gtype(gtype_t gtype){
    free(gtype);
    return;
}

/*ランダムにgtypeを作る*/
gtype_t mk_random_gtype(int code_length, int code_max){
    gtype_t ptr=mk_gtype(code_length);
    int i;
    for(i=0;i<code_length;i++){
        ptr[i]=rand()%(code_max+1);
    }
    return ptr;
}

/*実数とバイナリ・グレイ表現への変換*/

/*与えられた実数に近い、ビットストリング型のgtypeを作る*/
/*GRAY=1が指定されているときにはグレイ表現に変換*/
/*仕様：最後の行まで計算した後に残る端数は切り捨て*/
typedef int* gtype_t;
void encode_gtype(double value, gtype_t gtype,int code_length, double min, double max)
{double gap = max - min;
double remain_value = value-min;
            /*値のうち、遺伝子によって表現されている部分*/
double value_of_code;               /*その桁の遺伝子が表現する値*/
int position=1;
int pre_code=0;
int i=0;
int tmp;    /*グレイ表現変換用、一時保管変数*/
while(i<code_length){
    value_of_code=gap/pow(2,position);
    if(remain_value>=value_of_code){
        gtype[i]=1;
        remain_value -= value_of_code;
    }else{
        gtype[i] = 0;
    }
    /*グレイ表現への変換
     *バイナリ表現と、元のバイナリを右に１つシフトしたもののXORをとる
     */
    if(GRAY == 1){
        tmp = gtype[i];
        gtype[i] = (pre_code)^(gtype[i]);
        pre_code = tmp;
    }
    position++;
    i++;
}
return;
}

/*与えられたgtypeから実数に変換する方法*/
/*GRAY=1が指定されているときはグレイ表現に変換*/
double decode_gtype(gtype_t gtype, int code_length,
                    double min, double max)
{
    double gap = max - min;
    double decoded_value = min;
    int position = 1;
    int pre_code = 0;
        /*１つの上位の桁の表現(バイナリ),バイナリとグレイの変換に必要*/
    
    int i=0;
    /*グレイ表現の解釈*/
    /*変更されたバイナリの１つ上位の桁の表現との排他的論理和をとる*/
    if(GRAY==1){
        while (i<code_length){
            pre_code = pre_code^gtype[i];
            if(pre_code){
                decoded_value += gap/pow(2,position);
                /*最上位から順に、
                　最大値と最小値の差の1/2,1/4,1/8,1/16,,,となる*/
            }
            position++;
            i++;
        }
    }
    /*バイナリ表現の時*/
    else{
      while (i<code_length){
          if (gtype[i]){
              decoded_value += gap / pow(2,position);
              /*最上位から順に、
              　　最大値と最小値の差の1/2,1/4,1/8,1/16,,,となる*/
          }
          position++;
          i++;
      }
}
return decoded_value;
}

double decode_gtype2(gtype_t gtype, int code_length1, int code_length2,
                    double min, double max)
{
    double gap = max - min;
    double decoded_value = min;
    int position = 1;
    int pre_code = 0;
        /*１つの上位の桁の表現(バイナリ),バイナリとグレイの変換に必要*/
    
    int i=0;
    /*グレイ表現の解釈*/
    /*変更されたバイナリの１つ上位の桁の表現との排他的論理和をとる*/
    if(GRAY==1){
        while (i<code_length2){
            pre_code = pre_code^gtype[i+code_length1];
            if(pre_code){
                decoded_value += gap/pow(2,position);
                /*最上位から順に、
                　最大値と最小値の差の1/2,1/4,1/8,1/16,,,となる*/
            }
            position++;
            i++;
        }
    }
    /*バイナリ表現の時*/
    else{
      while (i<code_length2){
          if (gtype[i+code_length1]){
              decoded_value += gap / pow(2,position);
              /*最上位から順に、
              　　最大値と最小値の差の1/2,1/4,1/8,1/16,,,となる*/
          }
          position++;
          i++;
      }
}
return decoded_value;
}

void encode_gtype(double value, gtype_t gtype, int code_length, double min, double max);
double decode_gtype(gtype_t gtype, int code_length, double min, double max);

//

/*gtypeへのコピー*/
void copy_gtype(gtype_t new_gtype, gtype_t old_gtype, int length){
    int i=0;
    for(i=0;i<length;i++){
        new_gtype[i]=old_gtype[i];
    }
    return;
}

//int cross_gtype(gtype_t gtype1, gtype_t gtype2, int length);
//int mutate_gtype(gtype_t gtype, int length, int code_max, double pm);

/*交叉と突然変異の関数*/
typedef int* gtype_t;


/*-------------------1点交叉----------------------------------*/

/*
int cross_gtype(gtype_t gtype1, gtype_t gtype2, int length)
{
    int cross_point=rand()%(length-1);
    int i = cross_point+1;
    int tmp;
    while(i<length){
        tmp = gtype1[i];
        gtype1[i]=gtype2[i];
        gtype2[i]=tmp;
        i++;
    }
    return cross_point;
}
*/

/*-------------------一様交叉-----------------------------------*/

int cross_gtype(gtype_t gtype1, gtype_t gtype2, int length)
{
    int cross_point=0;
    //int i;
    int tmp;
    int k = 0;
    int r;

    while(k<length){
        r = rand()%(CODE_MAX+1);
        if(floor(r) == 1){
            tmp = gtype1[k];
            gtype1[k] = gtype2[k];
            gtype2[k] = tmp;
        }

        k++;
    }
    return cross_point;
}



/*-------------二点交叉----------------*/
/*
int cross_gtype(gtype_t gtype1, gtype_t gtype2, int length)
{
    int cross_point=rand()%(length-1);
    int cross_point2 = rand()%(length-1);
    int i = cross_point+1;
    int ii = cross_point2+1;
    int tmp;

    if (ii<i){
        while(ii<i){
            tmp = gtype1[ii];
            gtype1[ii] = gtype2[ii];
            gtype2[ii] = tmp;
            ii++;
        }
    }else if(i<ii){
        while(i<ii){
            tmp = gtype1[i];
            gtype1[i] = gtype2[i];
            gtype2[i] = tmp;
            i++;
        }

    }else{
    while(i<length){
        tmp = gtype1[i];
        gtype1[i]=gtype2[i];
        gtype2[i]=tmp;
        i++;
    }

    }

   return cross_point;
}
*/


int mutate_gtype(gtype_t gtype, int length, int code_max, double pm)
{
    //エラー処理
    if(pm>1.0 || pm<0.0){
        printf(
            "%f mutation probability must be from 0.0 to 1.0 \n",pm
        );
        exit(-1);
    }
    int mutate_point =0;
    int i=0;
    double rm;
    for(i=0;i<length;i++){
        rm=(double)rand()/RAND_MAX;
        if(rm<pm){
            gtype[i] = rand()% (code_max+1);
            mutate_point++;
        }
    }
    return mutate_point;
}

int cross_gtype(gtype_t gtype1, gtype_t gtype2, int length);
int mutate_gtype(gtype_t gtype, int lengthy, int code_max,double pm);

//

/*gtypeを表示する*/
void print_gtype(gtype_t gtype,int length){
    int i=0;
    printf("[");
    while(i<length){
        if(gtype[i]<10){
            printf("%d",gtype[i]);

        }else{
            printf("(%d)",gtype[i]);
        }
        i++;
    }
    printf("]");
}

/*線形リスト用の隣接した要素の入れ替え、引数は先頭のindivisu=al_tのアドレス*/
void switch_gene(indivisual_t *indivisual){
    indivisual_t tmp_ptr1=(*indivisual) ->next ->next;
    indivisual_t tmp_ptr2=(*indivisual) ->next;
    (*indivisual) ->next ->next =*indivisual;
    (*indivisual) ->next =tmp_ptr1;
    (*indivisual) =tmp_ptr2;
    return;
}

/*個体を作る。メモリ領域の確保。初期化*/
indivisual_t mk_gene(int code_length, int code_max){
    indivisual_t ptr = reinterpret_cast<indivisual_t>(my_malloc(sizeof(struct ga_indivisual)));
    ptr->gtype=mk_random_gtype(code_length, code_max);
    //ptr->ptype=0;
    ptr->ptype1=0;
    ptr->ptype2=0;
    ptr->ptype3=0;
    ptr->ptype4=0;
    ptr->ptype5=0;
    ptr->ptype6=0;
    //ptr->ptype7=0;
    ptr->fitness=0;
    ptr->next=NULL;
    ptr->parent1=0;
    ptr->parent2=0;
    ptr->cross_point=0;
    return ptr;
}

/*個体をコピーする*/
void copy_gene(indivisual_t new_gene, indivisual_t old_gene, int code_length){
    copy_gtype(new_gene->gtype,old_gene->gtype,code_length);
    new_gene->fitness=old_gene->fitness;
    new_gene->parent1=old_gene->rank;
    new_gene->parent2=old_gene->rank;
    new_gene->cross_point=code_length-1;
    return;
}

/*交叉、突然変異で子供を作る　突然変異を繰り返す*/
int mk_children_genes(indivisual_t child1, indivisual_t child2, indivisual_t parent1, indivisual_t parent2, int code_length,int code_max,double pm){
    int cross_point,mutateCount;
    copy_gene(child1,parent1,code_length);
    copy_gene(child2,parent2,code_length);
    cross_point=cross_gtype(child1->gtype,child2->gtype,code_length);
    child1->parent1=parent1->rank;
    child1->parent2=parent2->rank;
    child1->cross_point=cross_point;
    child2->parent1=parent2->rank;
    child2->parent2=parent1->rank;
    child2->cross_point=cross_point;
    mutateCount=mutate_gtype(child1->gtype,code_length,code_max,pm);
    mutateCount+=mutate_gtype(child2->gtype,code_length,code_max,pm);
    return mutateCount;
}

/*GA集団の作成、初期化を行う*/
ga_population_t mk_init_ga_population(int population_size, int code_length, int code_max){
    ga_population_t population=reinterpret_cast<ga_population_t>(my_malloc(sizeof(struct ga_population)));
    population->pselect=(double*)my_malloc(sizeof(double)*population_size);
    population->mutate_count=0;
    population->population_size=population_size;
    population->code_length=code_length;
    population->code_max=code_max;
    indivisual_t list_tale;
    population->genes=mk_gene(code_length, code_max);
    list_tale=population->genes;
    int i=1;
    for(i=1;i<population_size;i++){
        list_tale->next=mk_gene(code_length, code_max);
        list_tale=list_tale->next;
    }
    return population;

}

/*指定した文字chを指定した長さlengthだけ繰り返す関数*/
/*print_population(・)の中で使われる*/
void print_sequence(char ch, int length){
    int i=0;
    for(i=0;i<length;i++){
        printf("%c",ch);
}

}

/*集団を表示する*/
/*左から、世代数、親のインデックス、交叉点、gtype,ptype,fitnessを表示する*/
/*また、最後に突然変異の回数を表示する*/
void print_population(ga_population_t population){
    indivisual_t member=population->genes;
    int i=0;
    printf("---------");
    print_sequence('-',LENGTH+2);
    printf("-------\n");
    printf("# parents xsite gtype");
    print_sequence('-',LENGTH-3);
    printf("ptype fitness\n");//

    while(member != NULL){
        printf("%-3d (%3d,%3d) %3d ", i,member->parent1, member->parent2, member->cross_point);
        print_gtype(member->gtype,population->code_length);
        printf(" %+3.3f %3.3f   %3.3f   %3.3f   %3.3f   %3.3f   %3.3f\n",member->ptype1,member->ptype2,member->ptype3,member->ptype4,member->ptype5,member->ptype6,member->fitness);
        member=member->next;
        i++;
    }
    printf("total mutate %d\n",population->mutate_count);
    return;
}

/*適合度を出力
　最大、平均、最小
　CSV形式にする
*/
void print_fitness(ga_population_t population){
printf("%f,%f,%f,%f,%f,%f,%f,%f,%f",
    population->max_fitness,population->avg_fitness,
    population->min_fitness,population->genes->ptype1,population->genes->ptype2,population->genes->ptype3,population->genes->ptype4,population->genes->ptype5,population->genes->ptype6);
print_gtype(population->genes->gtype,
    population->code_length);
printf("\n");
return;
}

/*適合度を計算し、線形リストに適合度順に挿入する関数*/

/*適合度の比較関数*/
/*indivisualA　の適合度が小さければ１を返し、等しいか大きければ１を返す*/
int less_than(indivisual_t indivisualA, indivisual_t indivisualB)
{
    return (indivisualA->fitness<indivisualB->fitness);
}

/*適合度計算の関数*/
void calc_fitness(ga_population_t population,double value_min, double value_max, double m[])
{
    indivisual_t ptr = population->genes;
    indivisual_t next;
    indivisual_t indivisual_ptr = NULL;
    indivisual_t search_ptr = ptr;

    double x;
    int i=0;
    while(ptr != NULL){

        //適合度順に線形リストに挿入
        search_ptr = indivisual_ptr;
        if( search_ptr == NULL || less_than(indivisual_ptr, ptr)){
            ptr->next = indivisual_ptr;
            indivisual_ptr = ptr;
        }else{
            while(search_ptr->next != NULL){
                if(less_than(search_ptr->next, ptr)){
                    break;
                }
                search_ptr = search_ptr->next;
            }
            ptr->next = search_ptr->next;
            search_ptr->next = ptr;
        }
        ptr=next;
        i++;


    }
    population->genes = indivisual_ptr;

    return;
}

void calc_fitnessq(ga_population_t population)
{
    indivisual_t ptr = population->genes;
    indivisual_t next;
    indivisual_t indivisual_ptr = NULL;
    indivisual_t search_ptr = ptr;

    double x,y,z,k,d,a,b,c;
    //int i=0;
    while(ptr != NULL){
        x= decode_gtype(ptr->gtype,LENGTH1,MIN,MAX);
        y= decode_gtype2(ptr->gtype,LENGTH1,LENGTH2,MIN1,MAX1);
        z= decode_gtype2(ptr->gtype,LENGTH2,LENGTH3,MIN2,MAX2);
        a= decode_gtype2(ptr->gtype,LENGTH3,LENGTH4,MIN3,MAX3);
        b= decode_gtype2(ptr->gtype,LENGTH4,LENGTH5,MIN4,MAX4);
        c= decode_gtype2(ptr->gtype,LENGTH5,LENGTH6,MIN5,MAX5);

        ptr->ptype1 = x;//
        ptr->ptype2 = y;//
        ptr->ptype3 = z;
        ptr->ptype4 = a;//
        ptr->ptype5 = b;//
        ptr->ptype6 = c;
        k = F_X_Y_Z;//
        d=GG_K(k);
        ptr->fitness = d;
        next = ptr->next;
        ptr->next = NULL;

        search_ptr = indivisual_ptr;
        if( search_ptr == NULL || less_than(indivisual_ptr, ptr)){
            ptr->next = indivisual_ptr;
            indivisual_ptr = ptr;
        }else{
            while(search_ptr->next != NULL){
                if(less_than(search_ptr->next, ptr)){
                    break;
                }
                search_ptr = search_ptr->next;
            }
            ptr->next = search_ptr->next;
            search_ptr->next = ptr;
        }
        ptr=next;
        


    }
    population->genes = indivisual_ptr;
    
        //適合度順に線形リストに挿入

    return;
}


int less_than(indivisual_t indivisualA,indivisual_t indivisualB);
//void calc_fitness(ga_population_t population,double value_min, double value_max);

//

/*選択を実行する関数*/
void calc_pselect(ga_population_t population)
{
    int i;
    population->pselect[0] = population->genes->fitness;
    indivisual_t gene_ptr = population->genes->next;
    for (i=1; i < population->population_size; i++){
        population->pselect[i]=
            population->pselect[i-1] + gene_ptr->fitness;
        gene_ptr = gene_ptr->next;
    }
    for (i=0;i< population->population_size;i++){
        population->pselect[i] /=
            population->pselect[population->population_size-1];
    }
    return;
}

/*ルーレット方式による親選択*/
indivisual_t select_parent_roulette(ga_population_t population)
{
    int j=0;
    double r;
    indivisual_t parent;
    r = (double)rand()/RAND_MAX;
    parent = population->genes;
    while (r > population->pselect[j]){
        parent = parent->next;
        j++;
    }
    return parent;
}

/*トーナメント方式による親選択*/
indivisual_t select_parent_tournament(
    ga_population_t population, int tournament_size
)
{
    int pop = population->population_size;
    int i,j,r,min = pop;
    indivisual_t min_selected = NULL;
    indivisual_t ptr;
    for(i=0;i<tournament_size;i++){
        r = rand() % pop;
        if(min > r){
            min = r;
        }
    }
    ptr = population->genes;
    for(j=0;j<min;j++){
        ptr = ptr->next;
    }
    min_selected = ptr;
    return min_selected;
}

void calc_pselect(ga_population_t population);
indivisual_t select_parent_roulette(ga_population_t population);
indivisual_t select_parent_tournament(
    ga_population_t population, int tournament_size
);

//

/*親個体の選択、param.hのSELECTION_METHODによって
ルーレット選択かトーナメント選択を行う*/
indivisual_t select_parent(ga_population_t population){
    indivisual_t parent;
    switch(SELECTION_METHOD){
        case 1:
         parent=select_parent_roulette(population);
         break;
        case 2:
         parent=
            select_parent_tournament(population,TOURNAMENT_SIZE);
        break;
    default:
        fprintf(stderr,"invalid number on SELSCTION_METHOD\n");
        exit(1);
    }
    return parent;
}

/*適合度順に並んだ線形リストから
最大値、最小値、平均値を記録、順番付け*/
void normalize_population(ga_population_t population){
    int i;
    indivisual_t tmp;

    tmp=population->genes;
    population->max_fitness=population->genes->fitness;
                                /*先頭の適合度が最大適合度*/
    double avg=0.0;
    /*順番付け*/
    for(i=0;i<population->population_size;i++){
        avg += tmp->fitness;
        tmp->rank=i;
        if(tmp->next==NULL){
            population->avg_fitness=tmp->fitness;                                                                
        }/*最後尾の適合度が最小適合度*/
        tmp=tmp->next;
    }
    avg=avg/population->population_size;
    population->avg_fitness=avg;
    return;
}

/*新しい世代の生成
　new_populationのメモリ領域はすでに確保してあるとする
　必ずソート済みのpopulationを渡すこと
*/
void generate_population(ga_population_t new_population,ga_population_t old_population, double gap,double elete_rate, double mutate_prob,double crossover_prob)
{
    int num_of_remain=
(int)(old_population->population_size*(1-gap));
    /*親世代からコピーする数*/
    int num_of_elete=(int)(num_of_remain*elete_rate);
    /*コピー枠のうちのエリートの数*/
    int generated;
    double rand_double;
    indivisual_t old_gene=old_population->genes;
    indivisual_t new_gene=new_population->genes;
    /*親選択テーブルを準備*/
    calc_pselect(old_population);
    /*エリート戦略　親世代での上位一定数はそのまま子供になる*/
    for (generated=0;generated<num_of_elete;generated++){
        copy_gene(new_gene,old_gene,old_population->code_length);
        old_gene=old_gene->next;
        new_gene=new_gene->next;
    }
    /*エリート以外のそのまま子供になる枠*/
    for( ; generated<num_of_remain ; generated++){
        copy_gene(new_gene,select_parent(old_population),
                    old_population->code_length);
        new_gene=new_gene->next;
    }
    new_population->mutate_count=0;
    /*交叉・突然変異を適用する枠*/
    /*残りの個体数が奇数の時は、１つだけ突然変異で子供を作る*/
    if((old_population->population_size - generated)%2==1){
        copy_gene(new_gene,select_parent(old_population),
                  old_population->code_length);
        new_population->mutate_count+=
            mutate_gtype(new_gene->gtype,
                         old_population->code_length,
                         old_population->code_max,
                         mutate_prob);
        new_gene=new_gene->next;
        generated++;
    }
    /*交叉・突然をする*/
    for(; generated<old_population->population_size;generated+=2){
        rand_double=(double)rand()/RAND_MAX;
        /*交叉するとき*/
        if(rand_double<crossover_prob){
            new_population->mutate_count +=
                mk_children_genes(new_gene,
                new_gene->next,select_parent(old_population),
                select_parent(old_population),
                old_population->code_length,
                old_population->code_max,mutate_prob);
            new_gene=new_gene->next->next;
        }
        /*交叉しないとき*/
        else{
            copy_gene(new_gene,select_parent(old_population),
                      old_population->code_length);
            new_population->mutate_count+=
            mutate_gtype(new_gene->gtype,
            old_population->code_length,old_population->code_max,
            mutate_prob);
            new_gene=new_gene->next;
            copy_gene(new_gene,select_parent(old_population),
                      old_population->code_length);
            new_population->mutate_count+=
                mutate_gtype(new_gene->gtype,
                old_population->code_length,old_population->code_max,
                mutate_prob);
            new_gene=new_gene->next;
        }
    }
    return;
}

/*main関数*/
/*GAの実行*/
    int v=0;
    int r=0;
    int p=0;
int main(){
    double m[100000];
    double n[100000];
    double f[100000];
    //double et=2000;
    //double t_ini = 0;
    int v;
    int r;
    int p;

    /*乱数に引数を与える*/
    srand(time(NULL));
    ga_population_t parent_group=
        mk_init_ga_population(POP,LENGTH,CODE_MAX);
    ga_population_t child_group=
        mk_init_ga_population(POP,LENGTH,CODE_MAX);

    int i;
    if(PRINT_FITNESS==1){
        printf("#generation,max_fitness, avg_fitness,min_fitness, best_indivisual_ptype,best_indivisual_gtype\n");
    }
    for(i=0;i<=GENERATION;i++){
        //集団の適合度を計算し、線形リストを作る
        //calc_fitnessq(parent_group);
    
    indivisual_t ptr = parent_group->genes;
    indivisual_t next;
    indivisual_t indivisual_ptr = NULL;
    indivisual_t search_ptr = ptr;

    double x,y,z,k,d,a,b,c;

    while(ptr != NULL){
        x= decode_gtype(ptr->gtype,LENGTH1,MIN1,MAX1);
        y= decode_gtype2(ptr->gtype,LENGTH1,LENGTH2,MIN2,MAX2);
        z= decode_gtype2(ptr->gtype,LENGTH2,LENGTH3,MIN3,MAX3);
        a= decode_gtype2(ptr->gtype,LENGTH3,LENGTH4,MIN4,MAX4);
        b= decode_gtype2(ptr->gtype,LENGTH4,LENGTH5,MIN5,MAX5);
        c= decode_gtype2(ptr->gtype,LENGTH5,LENGTH6,MIN6,MAX6);
        ptr->ptype1 = x;
        ptr->ptype2 = y;
        ptr->ptype3 = z;
        ptr->ptype4 = a;
        ptr->ptype5 = b;
        ptr->ptype6 = c;

        //printf("%lf,%lf,%lf\n",x,y,z);
        k = F_X_Y_Z;//
        d=GG_K(k);

        //m[v]=d;
        /*n[r]=d;
        f[p]=d;*/

        /*printf("%lf\n",d);
        printf("%lf\n",m[v]);
        printf("%lf\n",n[r]);
        printf("%lf\n",f[p]);*/
 

        ptr->fitness = d;
        next = ptr->next;
        ptr->next = NULL;

        search_ptr = indivisual_ptr;
        if( search_ptr == NULL || less_than(indivisual_ptr, ptr)){
            ptr->next = indivisual_ptr;
            indivisual_ptr = ptr;
        }else{
            while(search_ptr->next != NULL){
                if(less_than(search_ptr->next, ptr)){
                    break;
                }
                search_ptr = search_ptr->next;
            }
            ptr->next = search_ptr->next;
            search_ptr->next = ptr;
        }
        //printf("%lf\n",m[v]);
        ptr=next;

        //v++;
        /*r++;
        p++;*/


    }
    parent_group->genes = indivisual_ptr;
    //printf("---------------------\n");





       

    //printf("---------------------\n");





        //最大値、最小値
        normalize_population(parent_group);

        //現在世代の表示
        if(PRINT_GROUP ==1){
            print_population(parent_group);

        }
        if(PRINT_FITNESS ==1){
            printf("%3d,",i);
            print_fitness(parent_group);

        }
        //現在世代parent_groupから次世代child_groupを作る
        generate_population(child_group,parent_group,GAP,
                            ELETE_RATE,P_MUTATE,P_CROSS);

        parent_group=child_group;

    }
    /*printf("----------------\n");
    for(v=0;v<=19;v++){
        printf("%lf\n",m[v]);
    }
    printf("----------------\n");
    for(r=0;r<=19;r++){
        printf("%lf\n",n[r]);
    }
    printf("----------------\n");
    for(p=0;p<=19;p++){
        printf("%lf\n",f[p]);
    }*/


    return 0;
}

