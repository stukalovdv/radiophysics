#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

class fdtd
{
    public:
        void displayLoading( int I , int N)           // Статус загрузки
        {
            cout << "Загрузка... " << ceil( I * 100 / N ) << "/" << N << "\% \r";
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


    private:
        const double c = 3e+10;             // Скорость света в СГС
        const double OMEGA_P_0 = 3e+9;      // Плазменная частота в вакууме до обезразмеривания
        //const double
        double NU_TILDA_;                   // Частота соударений
        double R2_TILDA_;                   // Радиус цилиндра
        double DELTA_;                      // Параметр неоднородности цилиндра
        double THETA_MULTIPLICATOR_;        // Множитель перед ПИ в угле наклона "чего-то там" (от 0 до 1/2)
        double KPD;                         // Эффективность излучения (КПД)

};

double simulation( double THETA_MULTIPLICATOR ,double NU_TILDA, double R2_TILDA, double DELTA ) // Функция вычислений
{
    double KPD = 0.871;

    int SIZE = 15;
    //vector <int> MASSIV(SIZE);
    vector <double> Er ( SIZE );
    vector <double> Ephi ( SIZE );
    vector <double> Ez ( SIZE );
    vector <double> Hr ( SIZE );
    vector <double> Hphi ( SIZE );
    vector <double> Hz ( SIZE );
    vector <double> Jr ( SIZE );
    vector <double> Jphi ( SIZE );
    vector <double> Jz ( SIZE );

    for ( int i = 0; i < SIZE; i++ )        // Начальные условия для полей E и H
    {
        Er[i] = 0;
        Ephi[i] = 0;
        Ez[i] = 0;
        Hr[i] = 0;
        Hphi[i] = 0;
        Hz[i] = 0;
    }
    return Er.size();
}




int main()
{
    setlocale(LC_ALL, "Russian");
    double THETA, NU_TILDA, R2_TILDA, DELTA;
    fdtd myCom;

    cout << "Запуск.\n";                        // Уведомление о запуске
    vector <int> a(10), b(10);
    cout << "Theta miltiplicator (| | x PI ) = ";  //
    cin >> THETA;                               //
    myCom.setTheta( THETA * 3.1415 );           //
    cout << "nu_tilda = ";                      //
    cin >> NU_TILDA;                            //        Ввод
    myCom.setNu( NU_TILDA );                    //
    cout << "R2_tilda = ";                      //      основных
    cin >> R2_TILDA;                            //
    myCom.setR2( R2_TILDA );                    //     параметров
    cout << "Delta = ";                         //
    cin >> DELTA;                               //       задачи
    myCom.setDelta( DELTA );                    //



    cout << simulation( myCom.getTheta() ,myCom.getNu(), myCom.getR2(), myCom.getDelta() ) << '\n';



    for (int i = 0; i < 100; i++)
    {
        myCom.displayLoading(i, 100);
    }


    myCom.displayFinish();                      // Уведомление о завершении программы
    cout << endl;
}
