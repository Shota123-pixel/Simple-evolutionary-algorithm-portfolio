typedef int* gtype_t;

typedef struct ga_indivisual* indivisual_t;
struct ga_indivisual{
    gtype_t gtype;//遺伝子型int配列の先頭ポインタ

    double ptype1;//表現型
    double ptype2;
    double ptype3;
    double ptype4;
    double ptype5;
    double ptype6;
    double fitness;//適合度

    indivisual_t next;//線形リストでの次の個体

    int rank;//線形リストの中での順位(ソート後)

    int parent1;//交叉での親1のインデックス
    int parent2;//交叉での親2のインデックス
    int cross_point;//交叉したポイント
};


typedef struct ga_population* ga_population_t;
struct ga_population{
    indivisual_t genes;//個体の線形リスト先頭へのポインタ
    double *pselect;//適合度の配列
    int mutate_count;//突然変異回数の合計
    double max_fitness;//適合度の最大値
    double min_fitness;//適合度の最小値
    double avg_fitness;//適合度の平均値
    int population_size;//集団の個体数
    int code_length;//遺伝子長
    int code_max;//各遺伝子座の最大値,ビットストリングの場合は1
};
/*
typedef int* gtype1_t;

typedef struct ga1_indivisual* indivisual1_t;
struct ga1_indivisual{
    gtype1_t gtype1;//遺伝子型int配列の先頭ポインタ
    double ptype1;//表現型
    double fitness1;//適合度
    indivisual1_t next1;//線形リストでの次の個体
    int rank1;//線形リストの中での順位(ソート後)
    int parent11;//交叉での親1のインデックス
    int parent21;//交叉での親2のインデックス
    int cross_point1;//交叉したポイント
};


typedef struct ga1_population* ga1_population_t;
struct ga1_population{
    indivisual1_t genes1;//個体の線形リスト先頭へのポインタ
    double *pselect1;//適合度の配列
    int mutate_count1;//突然変異回数の合計
    double max_fitness1;//適合度の最大値
    double min_fitness1;//適合度の最小値
    double avg_fitness1;//適合度の平均値
    int population_size1;//集団の個体数
    int code_length1;//遺伝子長
    int code_max1;//各遺伝子座の最大値,ビットストリングの場合は1
};*/
