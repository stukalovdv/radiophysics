#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

double fdtd( double THETA, double NU_TILDA, double R2_TILDA, double DELTA )
{
    double R1_tilda = R2_TILDA * ( 1 - DELTA );                     // Внутренний (обезразмеренный) радиус цилиндра
    double c = 3e+10, OMEGA_P_0 = 3e+9, J0 = 1;
    double NU = NU_TILDA * OMEGA_P_0;                               // Частота соударений 1

    double R1 = R1_tilda * c / OMEGA_P_0;                           // Внутренний радиус цилиндра
    double R2 = R2_TILDA * c / OMEGA_P_0;                           // Внешний радиус цилиндра
    double dr = 0.01 * NU_TILDA * ( R2 - R1 );                      // Шаг по пространству
    double dt = dr / ( c * 2 );                                     // Шаг по времени
    double T_MAX = 10;                                              // Расчетное (обезразмеренное) время
    int N_TIME = T_MAX / ( dt * OMEGA_P_0 );                        // Кол-во шагов по времени

    // Новые коэффициенты
    double V = c / sqrt( cos( THETA ) );                            //
    double MAIN_COEFFICIENT = c * dt / ( sin( THETA * THETA ) );    //
    double SUB_COEFFICIENT = c / V;                                 //

    //N_TIME = 5;
    vector <long double> T( N_TIME );                               // Вектор времени
    for ( int i = 0; i < N_TIME; i++ )                              //
    {                                                               //
        T[i] = dt * i;                                              //
    }                                                               //


    double R_MAX = R2 * 60;                                         // Расчетное пространство (по r)
    int NR = ( R_MAX + dr ) / dr;                                   // Кол-во шагов по пространству
    vector <double> r( NR ), r_tilda( NR ), r_alt( NR );            // Обозначаем r и 1/r
    for ( int i = 0; i < NR; i++ )                                  //
    {                                                               //
        r[i] = dr + dr * i;                                         //
        r_alt[i] = 1 / r[i];                                        //
        r_tilda[i] = r[i] * OMEGA_P_0 / c;                          //
    }                                                               //

    int NR1 =  R1 / dr, NR2 = R2 / dr;                              // Начало и конец неоднородного слоя (в кол-ве узлов)

    // Точка проверки значения поля от времени
    int FIELD_CHECK_POINT = NR2, FIELD_CHECK_POINT_TILDA = FIELD_CHECK_POINT * dr * OMEGA_P_0 / c;
    if ( NR1 == 0 ) NR1 = 1;


    vector <double> Fr( NR );                                       // Обозначаем F(r)
    for ( int i = 0; i < NR1; i++ )
    {
        Fr[i] = 1;
    }
    for ( int i = NR1; i < NR2; i++ )
    {
        Fr[i] = ( cos( (double)M_PI * (double)( i - NR1 ) / ( 2 * (int)( NR2 - NR1 ) ) ) ) * ( cos( (double)M_PI * (double)( i - NR1 ) / ( 2 * ( NR2 - NR1 ) ) ) );
    }
    for ( int i = NR2; i < NR; i++ )
    {
        Fr[i] = 0;
    }


    int N_PML = 10 * NR2;                                           // Поглощающий слой
    vector <double> sigma( NR );                                    // Начало поглощения
    for ( int i = 0; i < NR - N_PML; i++ )
    {
        sigma[i] = 1;
    }
    for ( int i = NR - N_PML; i < NR; i++ )
    {
        sigma[i] = cos( (double)( M_PI / 2 ) * ( i - ( NR - N_PML ) ) / ( NR - ( NR - N_PML ) ) );
        //sigma[i] = 1;
    }
    sigma[NR - 1] = 0;

    vector <double> Er( NR ), Ephi( NR ), Ez( NR );                     // Проекции электрических полей
    vector <double> Hr( NR ), Hphi( NR ), Hz( NR );                     // Проекции магнитных полей
    vector <double> Jr( NR ), Jphi( NR ), Jz( NR );                     // Проекции токовых компонент
    vector <double> Ept( N_TIME ), Hzt( N_TIME );                       // E_phi(t), H_z(t)

    //Начальные условия
    for ( int i = 0; i < NR; i++ )
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

    ofstream fout;
    string PATH = "main_data.dat";
    fout.open( PATH );
    fout << left << setw( 11 ) << "T" << "\t";
    fout << left << setw( 11 ) << "Er(t)" << "\t" << left << setw( 11 ) << "Ephi(t)" << "\t" << left << setw( 11 ) << "Ez(t)" << "\t";
    fout << left << setw( 11 ) << "Hr(t)" << "\t" << left << setw( 11 ) << "Hphi(t)" << "\t" << left << setw( 11 ) << "Hz(t)" << "\t";
    fout << left << setw( 11 ) << "Jr(t)" << "\t" << left << setw( 11 ) << "Jphi(t)" << "\t" << left << setw( 11 ) << "Jz(t)";
    fout << endl;

    for ( int n = 0; n < N_TIME; n++ )
    {
        //Jr
        for ( int i = 0; i < NR; i++ )
        {
            Jr[i] = Jr[i] * ( 1 - dt * NU ) + ( dt * OMEGA_P_0 * OMEGA_P_0 / ( 4 * M_PI ) ) * Fr[i] * Er[i];
        }
        //Jphi
        for ( int i = 0; i < NR; i++ )
        {
            Jphi[i] = Jphi[i] * ( 1 - dt * NU ) + ( dt * OMEGA_P_0 * OMEGA_P_0 / ( 4 * M_PI ) ) * Fr[i] * Ephi[i];
        }
        //Jz
        for ( int i = 0; i < NR; i++ )
        {
            Jz[i] = Jz[i] * ( 1 - dt * NU ) + ( dt * OMEGA_P_0 * OMEGA_P_0 / ( 4 * M_PI ) ) * Fr[i] * Ez[i];
        }
        //Er
        Er[0] = Er[1];
        for ( int i = 1; i < NR; i++ )
        {
            Er[i] = sigma[i] * ( Er[i] + MAIN_COEFFICIENT * ( r_alt[i] * Hz[i] + SUB_COEFFICIENT * ( Ez[i] - Ez[i - 1] ) / dr - ( 4 * M_PI / c ) * Jr[i] ) );
        }
        //Ephi
        Ephi[0] = Ephi[1];
        for ( int i = 1; i < NR; i++ )
        {
            Ephi[i] = sigma[i] * ( Ephi[i] + MAIN_COEFFICIENT * ( - SUB_COEFFICIENT * r_alt[i] * Ez[i] - ( Hz[i] - Hz[i - 1] ) / dr - ( 4 * M_PI / c ) * Jphi[i] ) );
        }
        //Ez
        Ez[0] = Ez[1];
        for ( int i = 1; i < NR; i++ )
        {
            Ez[i] = sigma[i] * ( Ez[i] + ( c * dt * r_alt[i] ) * ( ( Hphi[i] * r[i] - Hphi[i - 1] * r[i - 1] ) / dr - Hr[i] ) - 4 * M_PI * dt * Jz[i] );
        }
        Ez[NR - 1] = sigma[NR - 1] * ( Ez[NR - 1] - ( c * dt * r_alt[NR - 1] ) * Hr[NR - 1] );

        //Hr
        Hr[0] = Hr[1];
        for ( int i = 1; i < NR; i++ )
        {
            Hr[i] = sigma[i] * ( Hr[i] + MAIN_COEFFICIENT * ( r_alt[i] * Ez[i] + SUB_COEFFICIENT * ( ( Hz[i] - Hz[i - 1] ) / dr + ( 4 * M_PI / c ) * Jphi[i] ) ));
        }
        //Hphi
        Hphi[0] = Hphi[1];
        for ( int i = 1; i < NR; i++ )
        {
            Hphi[i] = sigma[i] * ( Hphi[i] + MAIN_COEFFICIENT * ( SUB_COEFFICIENT * ( r_alt[i] * Hz[i] - ( 4 * M_PI / c ) * Jr[i] ) + ( Ez[i] - Ez[i - 1] ) / dr ) );
        }
        //Hz
        Hz[0] = Hz[1];
        for ( int i = 1; i < NR; i++ )
        {
            Hz[i] = sigma[i] * ( Hz[i] - ( c * dt * r_alt[i] ) * ( Er[i] + ( r[i] * Ephi[i] - r[i - 1] * Ephi[i - 1] ) / dr ) );
        }

        Ept[n] = Ephi[FIELD_CHECK_POINT];
        Hzt[n] = Hz[FIELD_CHECK_POINT];

        cout << "Loading... " << ( n * 100 / N_TIME ) + 1 << "/" << 100 << "%\r";
        fout << left << setw( 11 ) << T[n] * OMEGA_P_0 << "\t";
        fout << left << setw( 11 ) << Er[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Ephi[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Ez[FIELD_CHECK_POINT] << "\t";
        fout << left << setw( 11 ) << Hr[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Hphi[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Hz[FIELD_CHECK_POINT] << "\t";
        fout << left << setw( 11 ) << Jr[( NR2 + NR1 ) / 2] << "\t" << left << setw( 11 ) << Jphi[( NR2 + NR1 ) / 2] << "\t" << left << setw( 11 ) << Jz[( NR2 + NR1 ) / 2] << "\t";
        fout << endl;
    }

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
    setlocale( LC_ALL,"Rus" ); // Русский язык в консоли

    // Основные параметры задачи
    double THETA_MULTIPLICATOR, NU_TILDA, R2_TILDA, DELTA;
    do
    {
        cout << "Theta miltiplicator (PI x | |) = ";
        cin >> THETA_MULTIPLICATOR;
    }
    while ( THETA_MULTIPLICATOR < 0 || THETA_MULTIPLICATOR > 0.5);
    cout << "nu = ";
    cin >> NU_TILDA;
    cout << "R2 = ";
    cin >> R2_TILDA;
    do
    {
        cout << "Delta = ";
        cin >> DELTA;
    }
    while ( DELTA < 0 || DELTA > 1 );

    cout << fdtd( THETA_MULTIPLICATOR * 3.1415926535, NU_TILDA, R2_TILDA, DELTA );
    cout << "\a\a\a\a\a\a\a\a\a\a\a\a\a\a" <<endl;


    //fstream
    //ifstream
    //ofstream


    //if (Fr.size() == r.size( ) ) cout << "\nYES";

    //cout << r.size() << "\t" << Fr.size() << "\t" << Nr << endl;




    return 0;
}
