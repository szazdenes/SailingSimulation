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
//    int roundedElevation = qRound(elevation);
    QFile file;

    if(ui->solRadioButton->isChecked()){
        if(ui->calciteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("./cal_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("./cal_sol_pm_ave.csv");
        }
        if(ui->cordieriteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("./cord_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("./cord_sol_pm_ave.csv");
        }
        if(ui->tourmalineCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("./tour_sol_am_ave.csv");
            if(time > 12)
                file.setFileName("./tour_sol_pm_ave.csv");
        }
    }

    if(ui->equRadioButton->isChecked()){
        if(ui->calciteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("./cal_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("./cal_equ_pm_ave.csv");
        }

        if(ui->cordieriteCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("./cord_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("./cord_equ_pm_ave.csv");
        }

        if(ui->tourmalineCheckBox->isChecked()){
            if(time <= 12)
                file.setFileName("./tour_equ_am_ave.csv");
            if(time > 12)
                file.setFileName("./tour_equ_pm_ave.csv");
        }

    }
    return northError;
}
