#include "sailingform.h"
#include "ui_sailingform.h"

SailingForm::SailingForm(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SailingForm)
{
    ui->setupUi(this);
}

SailingForm::~SailingForm()
{
    delete ui;
}

QVector2D SailingForm::getUnitStepVector(double Nerror)
{
    QVector2D result;
    result.setX(ui->speedDoubleSpinBox->value() * qCos(Nerror * M_PI / 180.0));
    result.setY(ui->speedDoubleSpinBox->value() * qSin(Nerror * M_PI / 180.0));
    return result;
}

double SailingForm::getNorthError(int time, int okta)
{
    double northError;
    QMap<int, double> elevationMap, roundedElevationMap;
//    int roundedElevation = qRound(elevation);
    QFile file;

    if(ui->solRadioButton->isChecked()){
        if(ui->calciteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("../cal_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("../cal_sol_pm_ave.csv");
        }
        if(ui->cordieriteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("../cord_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("../cord_sol_pm_ave.csv");
        }
        if(ui->tourmalineCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("../tour_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("../tour_sol_pm_ave.csv");
        }
        elevationMap = getTimeElevationMap("../elevation_Bergen_sol.dat");
    }

    if(ui->equRadioButton->isChecked()){
        if(ui->calciteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("../cal_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("../cal_equ_pm_ave.csv");
        }

        if(ui->cordieriteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("../cord_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("../cord_equ_pm_ave.csv");
        }

        if(ui->tourmalineCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("../tour_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("../tour_equ_pm_ave.csv");
        }
        elevationMap = getTimeElevationMap("../elevation_Bergen_equ.dat");
    }

    foreach(int key, elevationMap.keys()){
        if(ui->solRadioButton->isChecked() && qRound(elevationMap[key]) >= 50)
            roundedElevationMap[key] = 50.0;
        else if(ui->equRadioButton->isChecked() && qRound(elevationMap[key]) >= 25)
            roundedElevationMap[key] = 25.0;
        else
            roundedElevationMap[key] = qRound(elevationMap[key]);
    }

    QMap<QPair<int, int>, double> NErrorMap; /*first:elevation, second:okta*/

    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug("Opening error.");
    }
    QTextStream stream(&file);
    stream.readLine();
    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        int elev, cloud;
        double NError;
        linestream >> elev >> cloud >> NError;
        QPair<int, int> keyPair;
        keyPair.first = elev;
        keyPair.second = cloud;
        NErrorMap[keyPair] = NError;
    }

    return northError;
}

QMap<int, double> SailingForm::getTimeElevationMap(QString filename)
{
    QFile file(filename);
    QMap<int, double> timeElevMap;
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug("Opening error.");
    }
    QTextStream stream(&file);
    stream.readLine();
    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        int time;
        double elev;
        linestream >> time >> elev;
        timeElevMap[time] = elev;
    }
    return timeElevMap;
}

void SailingForm::on_startPushButton_clicked()
{
    getNorthError(5,2);
}
