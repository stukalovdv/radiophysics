// Новая задача
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

float c = 3e+10, OMEGA_P_0 = 3e+9, J0 = 1;                         // Глобальные переменные (скорость света, максимальная плазменная частота и плотность тока)

float fdtd( float THETA, float NU_TILDA, float R2_TILDA, float DELTA )
{
    float R1_tilda = R2_TILDA * ( 1 - DELTA );                     // Внутренний (обезразмеренный) радиус цилиндра
    float NU = NU_TILDA * OMEGA_P_0;                               // Частота соударений
    float R1 = R1_tilda * c / OMEGA_P_0;                           // Внутренний радиус цилиндра
    float R2 = R2_TILDA * c / OMEGA_P_0;                           // Внешний радиус цилиндра
    float dr = 0.001 * NU_TILDA * ( R2 - R1 );                     // Шаг по пространству
    float dt = dr / ( c * 20 );                                    // Шаг по времени
    float T_MAX = 2;                                               // Расчетное (обезразмеренное) время
    int N_TIME = T_MAX / ( dt * OMEGA_P_0 );                       // Кол-во шагов по времени
    //N_TIME = 25000;

    // Новые коэффициенты
    float MAIN_COEFFICIENT = c * dt / ( sin( THETA ) * sin( THETA ) ); // ПЕРЕД dt НУЖНО СТАВИТЬ с, НО ТОГДА БУДЕТ СЛИШКОМ СИЛЬНОЕ ОТКЛОНЕНИЕ И ОШИБКА
    float SUB_COEFFICIENT = cos( THETA );
    //float SUB_COEFFICIENT = 1.37e-35;
    if ( THETA > M_PI / 2 - 0.005 && THETA < M_PI / 2 + 0.005 )
    {
        SUB_COEFFICIENT = 0;
        MAIN_COEFFICIENT = c * dt;
    }

    //N_TIME = 25000;
    vector <float> T( N_TIME );                                     // Вектор времени
    for ( int i = 0; i < N_TIME; i++ )                              //
    {                                                               //
        T[i] = dt * i;                                              //
    }                                                               //


    float R_MAX = R2 * 40;                                         // Расчетное пространство (по r)
    int NR = ( R_MAX + dr ) / dr;                                   // Кол-во шагов по пространству
    vector <float> r( NR ), r_tilda( NR ), r_alt( NR );            // Обозначаем r и 1/r
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


    vector <float> Fr( NR );                                       // Обозначаем F(r)
    for ( int i = 0; i < NR1; i++ )
    {
        Fr[i] = 1;
    }
    for ( int i = NR1; i < NR2; i++ )
    {
        Fr[i] = ( cos( (float)M_PI * (float)( i - NR1 ) / ( 2 * (int)( NR2 - NR1 ) ) ) ) * ( cos( (float)M_PI * (float)( i - NR1 ) / ( 2 * ( NR2 - NR1 ) ) ) );
    }
    for ( int i = NR2; i < NR; i++ )
    {
        Fr[i] = 0;
    }


    int N_PML = 10 * NR2;                                           // Поглощающий слой
    vector <float> sigma( NR );                                    // Начало поглощения
    for ( int i = 0; i < NR - N_PML; i++ )
    {
        sigma[i] = 1;
    }
    for ( int i = NR - N_PML; i < NR; i++ )
    {
        sigma[i] = cos( (float)( M_PI / 2 ) * ( i - ( NR - N_PML ) ) / ( NR - ( NR - N_PML ) ) );
        //sigma[i] = 1;
    }
    sigma[NR - 1] = 0;

    vector <float> Er( NR ), Ephi( NR ), Ez( NR );                     // Проекции электрических полей
    vector <float> Hr( NR ), Hphi( NR ), Hz( NR );                     // Проекции магнитных полей
    vector <float> Jr( NR ), Jphi( NR ), Jz( NR );                     // Проекции токовых компонент
    vector <float> Ept( N_TIME ), Hzt( N_TIME );                       // E_phi(t), H_z(t)
    vector <float> Ezt( N_TIME ), Hpt( N_TIME );                       // E_z(t), H_phi(t)
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
        Jphi[i] = - J0 * Fr[i];
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
        /*
        //Jr
        for ( int i = 0; i < NR; i++ )
        {
            Jr[i] = J0 * sin( OMEGA_P_0 * n * dt );
        }
        //Jphi
        for ( int i = 0; i < NR; i++ )
        {
            Jphi[i] = - J0 * sin( OMEGA_P_0 * n * dt );
        }*/

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
            Er[i] = sigma[i] * ( Er[i] + MAIN_COEFFICIENT * ( r_alt[i] * Hz[i] - ( 4 * M_PI / c ) * Jr[i] + SUB_COEFFICIENT * ( Ez[i] - Ez[i - 1] ) / dr ) );
        }
        //Ephi
        Ephi[0] = Ephi[1];
        for ( int i = 1; i < NR; i++ )
        {
            Ephi[i] = sigma[i] * ( Ephi[i] - ( c * dt / dr ) * ( Hz[i] - Hz[i-1] ) - ( 4 * M_PI * dt ) * Jphi[i] );
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
            Hr[i] = sigma[i] * ( Hr[i] + MAIN_COEFFICIENT * ( r_alt[i] * Ez[i] + SUB_COEFFICIENT * ( ( Hz[i] - Hz[i - 1] ) / dr + ( 4 * M_PI / c ) * Jphi[i] ) ) );
        }
        //Hphi
        Hphi[0] = Hphi[1];
        for ( int i = 1; i < NR; i++ )
        {
            Hphi[i] = sigma[i] * ( Hphi[i] + MAIN_COEFFICIENT * ( SUB_COEFFICIENT * ( r_alt[i] * Hz[i] - ( 4 * M_PI / c ) * Jr[i] ) + ( Ez[i] - Ez[i - 1] ) / dr ) );
        }
        //Hz
        for ( int i = 0; i < NR - 1; i++ )
        {
            Hz[i] = sigma[i] * ( Hz[i] - ( c * dt * r_alt[i] ) * ( Er[i] + ( r[i + 1] * Ephi[i + 1] - r[i] * Ephi[i] ) / dr ) );
        }
        Hz[ NR - 1 ] = sigma[ NR - 1 ] * ( Hz[ NR - 1 ] - ( c * dt * r_alt[ NR - 1 ] ) * Er[ NR - 1 ]);


        Ept[n] = Ephi[FIELD_CHECK_POINT];
        Hzt[n] = Hz[FIELD_CHECK_POINT];

        Ezt[n] = Ez[FIELD_CHECK_POINT];
        Hpt[n] = Hphi[FIELD_CHECK_POINT];

        cout << "Loading... " << ( n * 100 / N_TIME ) + 1 << "/" << 100 << "%\r";
        fout << left << setw( 11 ) << T[n] * OMEGA_P_0 << "\t";
        fout << left << setw( 11 ) << Er[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Ephi[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Ez[FIELD_CHECK_POINT] << "\t";
        fout << left << setw( 11 ) << Hr[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Hphi[FIELD_CHECK_POINT] << "\t" << left << setw( 11 ) << Hz[FIELD_CHECK_POINT] << "\t";
        fout << left << setw( 11 ) << Jr[( NR2 + NR1 ) / 2] << "\t" << left << setw( 11 ) << Jphi[( NR2 + NR1 ) / 2] << "\t" << left << setw( 11 ) << Jz[( NR2 + NR1 ) / 2] << "\t";
        fout << endl;
    }

    fout.close();
    float I = 0;
    I += ( Ept[0] * Hzt[0] + Ept[ N_TIME-1 ] * Hzt[ N_TIME - 1 ] ) / 2;
    for ( int i = 1; i < N_TIME; i++ )
    {
        I += ( Ept[i] * Hzt[i] );
    }
    I += ( Ezt[0] * Hpt[0] + Ezt[ N_TIME-1 ] * Hpt[ N_TIME - 1 ] ) / 2;
    for ( int i = 1; i < N_TIME; i++ )
    {
        I += ( Ezt[i] * Hpt[i] );
    }
    float W_zap, W_izl;
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
    float THETA_MULTIPLICATOR, NU_TILDA, R2_TILDA, DELTA;
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

    cout << fdtd( THETA_MULTIPLICATOR * M_PI, NU_TILDA, R2_TILDA, DELTA );
    cout << "\a" <<endl;

    //cout << THETA_MULTIPLICATOR * M_PI << "  " << cos( THETA_MULTIPLICATOR * M_PI );

    //fstream
    //ifstream
    //ofstream


    //if (Fr.size() == r.size( ) ) cout << "\nYES";

    //cout << r.size() << "\t" << Fr.size() << "\t" << Nr << endl;




    return 0;
}
