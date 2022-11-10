/*GAで使われるパラメータをまとめて記述するファイル*/

#define F_X_Y_Z x+y+z+a+b+c-6 //探索する関数。関数への入力はx,関数からの出力はy
#define G_K 1/(1+fabs(k)) /*関数の出力値に対する適合度。適合度は高ければ高いほど良いものと扱われる。f(x)=0が最良*/
#define GRAY 1 /*グレイ表現なら1,バイナリ表現なら0*/
/*GAのパラメータ群*/

#define MAX (2.00)/*扱う実数の最大値*/
#define MAX1 (2.00)
#define MAX2 (2.00)
#define MAX3 (2.00)
#define MAX4 (2.00)
#define MAX5 (2.00)
#define MAX6 (2.00)
#define MIN (-2.00)/*扱う実数の最小値*/
#define MIN1 (-2.00)
#define MIN2 (-2.00)
#define MIN3 (-2.00)
#define MIN4 (-2.00)
#define MIN5 (-2.00)
#define MIN6 (-2.00)

#define LENGTH (180)/*遺伝子のコード長*/
#define LENGTH1 (30)
#define LENGTH2 (60)
#define LENGTH3 (90)
#define LENGTH4 (120)
#define LENGTH5 (150)
#define LENGTH6 (180)

//#define LENGTH ()

#define POP 100/*個体数*/
//#define POPR 100
#define CODE_MAX 1/*各遺伝子コードの最大値。これが1ならコードは0か1になる。ビット文字列の場合は１*/
//#define CODE_MAX1 1
#define GAP 0.9/*1回の生殖で子供と入れ替わる割合*/
//#define GAP1 0.9
#define ELETE_RATE 1.0/*そのまま生き残る数のうち、エリートの割合*/
//#define ELETE_RATE1 1.0
#define P_MUTATE 0.0333/*突然変異率。LENGTHの逆数程度が良い*/
//#define P_MUTATE1 0.0333
#define P_CROSS 1.0/*交叉確率*/
//#define P_CROSS1 1.0
#define GENERATION 200/*GAを計算する世代数*/
#define SELECTION_METHOD 2/*１はルーレット、２はトーナメント*/

/*トーナメントサイズ。トーナメントの時だけ意味がある*/
#define TOURNAMENT_SIZE 2

/*出力*/
#define PRINT_GROUP 1
#define PRINT_FITNESS 1