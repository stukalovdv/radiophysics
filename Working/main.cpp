#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

const double c = 3e+10;
const double OMEGA_P_0 = 3e+9;
const double PI = 3.14159265358979;

class fdtd
{
    public:
        void displayLoading( int I , int N)  // Статус загрузки
        {
            cout << "Загрузка... " << I * 100 / ( N - 1 ) << "/" << 100 << "\% \r";
        }
        void displayFinish()            // Статус завершения
        {
            cout << "\nWell done!\a";
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
    private:
        double NU_TILDA_;                   // Частота соударений
        double R2_TILDA_;                   // Радиус цилиндра
        double DELTA_;                      // Параметр неоднородности цилиндра
        double THETA_MULTIPLICATOR_;        // Множитель перед ПИ в угле наклона "чего-то там" (от 0 до 1/2)
};


double simulation( double THETA_MULTIPLICATOR, double NU_TILDA, double R2_TILDA, double DELTA ) // Функция вычислений
{

    fdtd myCom;

    double R1_TILDA = R2_TILDA* ( 1 - DELTA );

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

    int FIELD_CHECK_POINT = NR2 + 1;

    // Задаем поглощающий слой
    vector <double> sigma( NR );
    for ( int i = 0; i < NR; i++ )
    {
        sigma[i] = 1;
    }


    //setConstants( THETA_MULTIPLICATOR, NU_TILDA, R2_TILDA, DELTA );
    //double KPD = 0.871;

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
        Hz[i] = 0;
        Jr[i] = J0 * Fr[i];
        Jphi[i] = -J0 * Fr[i];
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

    vector <double> Ept( N_TIME ), Hzt( N_TIME );
    ofstream fout;
    string PATH = "superfile.dat";
    fout.open( PATH );
    fout << left << setw( 11 ) << "T" << "\t";
    fout << left << setw( 11 ) << "Er(t)" << "\t" << left << setw( 11 ) << "Ephi(t)" << "\t" << left << setw( 11 ) << "Ez(t)" << "\t";
    fout << left << setw( 11 ) << "Hr(t)" << "\t" << left << setw( 11 ) << "Hphi(t)" << "\t" << left << setw( 11 ) << "Hz(t)" << "\t";
    fout << left << setw( 11 ) << "Jr(t)" << "\t" << left << setw( 11 ) << "Jphi(t)" << "\t" << left << setw( 11 ) << "Jz(t)";
    fout << endl;

    for (int n = 0; n < N_TIME - 1; n++){
        //Jr
        for (int i = 0; i < r.size(); i++)
        {
            Jr[i] = Jr[i] * (1 - dt * NU) + (dt * OMEGA_P_0 * OMEGA_P_0 / (4 * PI)) * Fr[i] * Er[i];
        }
        //Jp
        for (int i = 0; i < r.size(); i++)
        {
            Jphi[i] = Jphi[i] * (1 - dt * NU) + (dt * OMEGA_P_0 * OMEGA_P_0 / (4 * PI)) * Fr[i] * Ephi[i];
        }
        //Er
        Er[0] = Er[1];
        for (int i = 1; i < r.size(); i++)
        {
            Er[i] = sigma[i] * (Er[i] + (c * dt * r_alt[i]) * Hz[i] - (4 * PI * dt) * Jr[i]);
        }
        //Ephi
        Ephi[0] = Ephi[1];
        for (int i = 1; i < r.size(); i++)
        {
            Ephi[i] = sigma[i] * (Ephi[i] - (c * dt / dr) * (Hz[i] - Hz[i-1]) - (4 * PI * dt) * Jphi[i]);
        }
        //Hz

        for (int i = 0; i < r.size() - 1; i++)
        {
            Hz[i] = sigma[i] * (Hz[i] - (c * dt * r_alt[i]) * Er[i] - (c * dt * r_alt[i] / dr) * (r[i + 1] * Ephi[i + 1] - r[i] * Ephi[i]));
        }
        Hz[r.size() - 1] = sigma[r.size() - 1] * (Hz[r.size() - 1] - (c * dt * r_alt[r.size() - 1]) * Er[r.size() - 1]);

        Ept[n] = Ephi[FIELD_CHECK_POINT];
        Hzt[n] = Hz[FIELD_CHECK_POINT];

        cout << "Loading... " << ( n * 100 / N_TIME ) + 1 << "/" << 100 << "%\r";
        fout << left << setw( 11 ) << T[n] * OMEGA_P_0 << "\t";
        fout << left << setw( 11 ) << Er[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Ephi[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Ez[FIELD_CHECK_POINT] << "\t";
        fout << left << setw( 11 ) << Hr[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Hphi[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Hz[FIELD_CHECK_POINT] << "\t";
        fout << left << setw( 11 ) << Jr[( NR2 + NR1 ) / 2] << "\t" << left << setw( 11 ) << Jphi[( NR2 + NR1 ) / 2] << "\t" << left << setw( 11 ) << Jz[( NR2 + NR1 ) / 2] << "\t";
        fout << endl;
    }
    /*
    for ( int i = 0; i < NR; i++ )
    {
        fout << Hz[i] << endl;
        //fout << Er[i] << "\t" << Ephi[i] << "\t" << Ez[i] << "\t" << Jr[i] << "\t" << Jphi[i] << "\t" << Jz[i] << "\t" << Fr[i] <<endl;
    }
    */
    fout.close();
    double I = 0;
    for ( int i = 1; i < N_TIME; i++ )
    {
        I += ( Ept[i] * Hzt[i] );
    }
    I += ( Ept[0] * Hzt[0] + Ept[N_TIME-1] * Hzt[N_TIME - 1] ) / 2;
    double W_zap, W_izl;
    W_zap = ( R2 * R2 + R1 * R1 ) * ( ( M_PI * J0 ) / OMEGA_P_0 ) * ( ( M_PI * J0 ) / OMEGA_P_0 );
    W_izl = c * ( dt / 4 ) * FIELD_CHECK_POINT * dr * I;

    cout << "\rWell done!         \n";
    cout << "File saved in path: " << PATH << endl;
    return W_izl / W_zap;
}




int main()
{
    setlocale(LC_ALL, "Rus");
    double THETA = 1, NU_TILDA, R2_TILDA, DELTA = 0;
    fdtd myCom;

    cout << "Запуск.\n";                        // Уведомление о запуске
    vector <int> a(10), b(10);
    while ( THETA > 0.5 || THETA < 0)
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



    cout << "\n" << simulation( myCom.getTheta() ,myCom.getNu(), myCom.getR2(), myCom.getDelta() ) << '\n';



    for (int i = 0; i < 100; i++)
    {
        myCom.displayLoading(i, 100);
    }


    myCom.displayFinish();                      // Уведомление о завершении программы
    cout << endl;
}
