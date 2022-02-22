#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

class fdtd
{
    public:
        void displayLoading( int I , int N)  // Статус загрузки
        {
            cout << "Загрузка... " << I * 100 / ( N - 1 ) << "/" << N << "\% \r";
        }
        void displayFinish()            // Статус завершения
        {
            cout << "\nWell done!";
        }
        void setNu( double NU_TILDA )        // Присваивание "Nu" с чертой
        {
            NU_TILDA_ = NU_TILDA;
        }
        void setR2( double R2_TILDA )        // Присваивание "R2" с чертой
        {
            R2_TILDA_ = R2_TILDA;
        }
        void setDelta( double DELTA )       // Присваивание "Delta" с чертой
        {
            DELTA_ = DELTA;
        }
        void setTheta( double THETA )       // Присваивание "Theta"
        {
            THETA_MULTIPLICATOR_ = THETA;
        }

        double getTheta()                   // Возврат можителя перед ПИ в "Theta"
        {
            return THETA_MULTIPLICATOR_;
        }

        double getNu()                      // Возврат "Nu" с чертой
        {
            return NU_TILDA_;
        }
        double getR2()                      // Возврат "R2" с чертой
        {
            return R2_TILDA_;
        }
        double getDelta()                   // Возврат "Delta" с чертой
        {
            return DELTA_;
        }
        double getC()
        {
            return c;
        }
        double getOmega()
        {
            return OMEGA_P_0;
        }


    private:
        const double c = 3e+10;             // Скорость света в СГС
        const double OMEGA_P_0 = 3e+9;      // Плазменная частота в вакууме до обезразмеривания
        //const double
        double NU_TILDA_;                   // Частота соударений
        double R2_TILDA_;                   // Радиус цилиндра
        double DELTA_;                      // Параметр неоднородности цилиндра
        double THETA_MULTIPLICATOR_;        // Множитель перед ПИ в угле наклона "чего-то там" (от 0 до 1/2)
};


double simulation( double THETA_MULTIPLICATOR ,double NU_TILDA, double R2_TILDA, double DELTA ) // Функция вычислений
{
    fdtd myCom;

    double R1_TILDA = R2_TILDA* ( 1 - DELTA );

    const double c = 3e+10;
    const double OMEGA_P_0 = 3e+9;
    const double PI = 3.14159265358979;

    double NU = NU_TILDA * OMEGA_P_0, R1 = R1_TILDA * c / OMEGA_P_0, R2 = R2_TILDA * c / OMEGA_P_0;
    double J0 = 1;

    double dr = 0.01 * NU_TILDA * ( R2 - R1 );
    double dt = dr / ( 2 * c );

    double T_MAX = 20 / OMEGA_P_0;
    int N_TIME = T_MAX / dt;

    double R_MAX = R2 * 60;
    int NR = R_MAX / dr;

    int NR1 = ceil( R1 / dr );
    int NR2 = ceil( R2 / dr );

    int N_PML = 10 * NR2;



    // Определения F(r)

    vector <double> Fr( NR );
    for ( int i = 0; i < NR1; i++ )
    {
        Fr[i] = 1;
    }
    for ( int i = NR1; i < NR2; i++ )
    {
        Fr[i] = (cos((double) PI * (double)( i - NR1 ) / ( 2 * (int)( NR2 - NR1 ) ))) * (cos((double) PI * (double)( i - NR1 ) / ( 2 * ( NR2 - NR1 ) )));
    }
    for ( int i = NR2; i < NR; i++)
    {
        Fr[i] = 0;
    }

    int FIELD_CHECK_POINT = NR2;

    /* Задаем поглощающий слой





    */

    //setConstants( THETA_MULTIPLICATOR, NU_TILDA, R2_TILDA, DELTA );
    double KPD = 0.871;

    vector <double> Er ( NR );
    vector <double> Ephi ( NR );
    vector <double> Ez ( NR );
    vector <double> Hr ( NR );
    vector <double> Hphi ( NR );
    vector <double> Hz ( NR );
    vector <double> Jr ( NR );
    vector <double> Jphi ( NR );
    vector <double> Jz ( NR );
    vector <double> T ( N_TIME );
    for ( int i = 0; i < NR; i++ )        // Начальные условия
    {
        Er[i] = 0;
        Ephi[i] = 0;
        Ez[i] = 0;
        Hr[i] = 0;
        Hphi[i] = 0;
        Hz[i] = 0;
        Jr[i] = J0 * Fr[i];
        Jphi[i] = -J0 * Fr[i];
        Jz[i] = 0;
    }
    for ( int i = 0; i < N_TIME; i++ )
    {
        T[i] = i * dt;
    }
    vector <double> r( NR ), r_alt( NR );
    for ( int i = 0; i < NR; i++ )
    {
        r[i] = i * dr;
        r_alt[i] = 1 / ( i * dr );
    }

    vector <double> Ert( N_TIME ), Hrt( N_TIME ), Ept( N_TIME );

    /*
    ofstream fout;
    fout.open("myFile.dat");
    for (int i = 0; i < NR; i++)
    {
        fout << Fr[i] << endl;
        cout << Fr[i] << endl;
    }
    */





    return KPD;
}




int main()
{
    //setlocale(LC_ALL, "Russian");
    double THETA = 1, NU_TILDA, R2_TILDA, DELTA = 0;
    fdtd myCom;

    cout << "Запуск.\n";                        // Уведомление о запуске
    vector <int> a(10), b(10);
    while ( THETA > 0.5 )
    {
        cout << "Theta miltiplicator (| | x PI ) = ";
        cin >> THETA;                           //
    }
    myCom.setTheta( THETA * 3.1415 );           //
    cout << "nu_tilda = ";                      //
    cin >> NU_TILDA;                            //        Ввод
    myCom.setNu( NU_TILDA );                    //
    cout << "R2_tilda = ";                      //
    cin >> R2_TILDA;                            //      основных
    myCom.setR2( R2_TILDA );                    //
    while ( DELTA <= 0 || DELTA > 1)            //     параметров
    {
        cout << "Delta = ";                                                                   //
        cin >> DELTA;                           //
    }                                           //       задачи
    myCom.setDelta( DELTA );                    //



    cout << simulation( myCom.getTheta() ,myCom.getNu(), myCom.getR2(), myCom.getDelta() ) << '\n';



    for (int i = 0; i < 100; i++)
    {
        myCom.displayLoading(i, 100);
    }


    myCom.displayFinish();                      // Уведомление о завершении программы
    cout << endl;
}
