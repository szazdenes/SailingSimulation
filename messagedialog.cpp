#include "messagedialog.h"
#include "ui_messagedialog.h"

MessageDialog::MessageDialog(QString message, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MessageDialog)
{
    ui->setupUi(this);
    ui->label->setText(message);
}

MessageDialog::~MessageDialog()
{
    delete ui;
}
