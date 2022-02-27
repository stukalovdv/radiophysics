#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;


double fdtd( double NU_TILDA, double R2_TILDA, double DELTA )
{
    double R1_tilda = R2_TILDA * ( 1 - DELTA );                     // Внутренний (обезразмеренный) радиус цилиндра
    double PI = 3.1415; // Пи
    double c = 3e+10, OMEGA_P_0 = 3e+9, J0 = 1;
    double NU = NU_TILDA * OMEGA_P_0;                               // Частота соударени
    double R1 = R1_tilda * c / OMEGA_P_0;                           // Внутренний радиус цилиндра
    double R2 = R2_TILDA * c / OMEGA_P_0;                           // Внешний радиус цилиндра
    double dr = 0.01 * NU_TILDA * ( R2 - R1 );                      // Шаг по пространству
    double dt = dr / ( c * 2 );                                     // Шаг по времени
    double T_MAX = 40;                                              // Расчетное (обезразмеренное) время
    int N_TIME = T_MAX / ( dt * OMEGA_P_0 );                        // Кол-во шагов по времени


    // Вектор времени
    vector <double> T( N_TIME );
    for ( int i = 0; i < N_TIME; i++ )
    {
        T[i] = dt * i;
    }

    // Обозначаем r и 1/r
    double R_MAX = R2 * 60;                                         // Расчетное пространство (по r)
    int NR = ( R_MAX + dr ) / dr;                                   // Кол-во шагов по пространству
    vector <double> r( NR ), r_tilda( NR ), r_alt( NR );
    for ( int i = 0; i < NR; i++ )
    {
        r[i] = dr + dr * i;
        r_alt[i] = 1 / r[i];
        r_tilda[i] = r[i] * OMEGA_P_0 / c;
    }

    int NR1 =  R1 / dr, NR2 = R2 / dr , N_PML = 10 * NR2;
    int FIELD_CHECK_POINT = NR2, FIELD_CHECK_POINT_TILDA = FIELD_CHECK_POINT * dr * OMEGA_P_0 / c;
    if ( NR1 == 0 ) NR1 = 1;

    // Обозначаем F(r)
    vector <double> Fr( NR );
    for ( int i = 0; i < NR1; i++ )
    {
        Fr[i] = 1;
    }
    for ( int i = NR1; i < NR2; i++ )
    {
        Fr[i] = ( cos( (double )PI * (double)( i - NR1 ) / ( 2 * (int)( NR2 - NR1 ) ) ) ) * ( cos( (double)PI * (double)( i - NR1 ) / ( 2 * ( NR2 - NR1 ) ) ) );
    }
    for ( int i = NR2; i < NR; i++ )
    {
        Fr[i] = 0;
    }

    // PML
    vector <double> sigma( NR );
    for ( int i = 0; i < NR - N_PML; i++ )
    {
        sigma[i] = 1;
    }
    for ( int i = NR - N_PML; i < NR; i++ )
    {
        sigma[i] = cos( (double)( PI / 2 ) * ( i - ( NR - N_PML ) ) / ( NR - ( NR - N_PML ) ) );
        //sigma[i] = 1;
    }
    sigma[NR - 1] = 0;

    vector <double> Er( NR ), Ephi( NR ), Ez( NR );
    vector <double> Hr( NR ), Hphi( NR ), Hz( NR );
    vector <double> Jr( NR ), Jphi( NR ), Jz( NR );
    vector <double> Ept( N_TIME ), Hzt( N_TIME );

    //íà÷àëüíûå óñëîâèÿ
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
    fout.open( "../python/file.txt" );

    //âû÷èñëåíèå
    for ( int n = 0; n < N_TIME; n++ )
    {
        //Jr
        for ( int i = 0; i < NR; i++ )
        {
            Jr[i] = Jr[i] * ( 1 - dt * NU ) + ( dt * OMEGA_P_0 * OMEGA_P_0 / ( 4 * PI ) ) * Fr[i] * Er[i];
        }
        //Jphi
        for ( int i = 0; i < NR; i++ )
        {
            Jphi[i] = Jphi[i] * ( 1 - dt * NU ) + ( dt * OMEGA_P_0 * OMEGA_P_0 / ( 4 * PI ) ) * Fr[i] * Ephi[i];
        }
        //Jz
        for ( int i = 0; i < NR; i++ )
        {
            Jz[i] = Jz[i];
        }
        //Er
        Er[0] = Er[1];
        for ( int i = 1; i < NR; i++ )
        {
            Er[i] = sigma[i] * ( Er[i] + ( c * dt * r_alt[i] ) * Hz[i] - ( 4 * PI * dt ) * Jr[i] );
        }
        //Ephi
        Ephi[0] = Ephi[1];
        for ( int i = 1; i < NR; i++ )
        {
            Ephi[i] = sigma[i] * ( Ephi[i] - ( c * dt / dr ) * ( Hz[i] - Hz[i-1] ) - ( 4 * PI * dt ) * Jphi[i] );
        }
        //Ez
        for ( int i = 0; i < NR; i++ )
        {
            Ez[i] = Ez[i];
        }

        //Hr
        for ( int i = 0; i < NR; i++ )
        {
            Hr[i] = Hr[i];
        }
        //Hphi
        for ( int i = 0; i < NR; i++ )
        {
            Hphi[i] = Hphi[i];
        }
        //Hz
        for ( int i = 0; i < NR - 1; i++ )
        {
            Hz[i] = sigma[i] * ( Hz[i] - ( c * dt * r_alt[i] ) * Er[i] - ( c * dt * r_alt[i] / dr ) * ( r[i + 1] * Ephi[i + 1] - r[i] * Ephi[i] ) );
        }
        Hz[NR - 1] = sigma[NR - 1] * ( Hz[NR - 1] - ( c * dt * r_alt[NR - 1] ) * Er[NR - 1] );

        Ept[n] = Ephi[FIELD_CHECK_POINT];
        Hzt[n] = Hz[FIELD_CHECK_POINT];

        cout << "Loading... " << ( n * 100 / N_TIME ) + 1 << "/" << 100 << "%\r";

        fout << Hz[FIELD_CHECK_POINT] << endl;
    }

    double I = 0;
    for ( int i = 1; i < N_TIME; i++ )
    {
        I += ( Ept[i] * Hzt[i] );
    }
    I += ( Ept[0] * Hzt[0] + Ept[N_TIME-1] * Hzt[N_TIME - 1] ) / 2;
    double W_zap, W_izl;
    W_zap = ( R2 * R2 + R1 * R1 ) * ( ( PI * J0 ) / OMEGA_P_0 ) * ( ( PI * J0 ) / OMEGA_P_0 );
    W_izl = c * ( dt / 4 ) * FIELD_CHECK_POINT * dr * I;


    if ( !fout.is_open() )
    {
        cout << "Ошибка записи файла!" << endl;
    }
    /*
    else{
        for ( int i = 0; i < Hzt.size(); i++ )
        {
            //fout << setprecision(20) << fixed;
            fout << Hzt[i] << " ";
        }
    }
    */

    fout.close();
    cout << "\rWell done!         \n";

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

    cout << fdtd( NU_TILDA, R2_TILDA, DELTA );
    cout << "\a\a\a\a\a\a\a\a\a\a\a\a\a\a" <<endl;


    //fstream
    //ifstream
    //ofstream


    //if (Fr.size() == r.size( ) ) cout << "\nYES";

    //cout << r.size() << "\t" << Fr.size() << "\t" << Nr << endl;




    return 0;
}
