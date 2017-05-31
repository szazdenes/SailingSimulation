#include "contourrecognition.h"

ContourRecognition::ContourRecognition(QObject *parent) : QObject(parent)
{

}

ContourRecognition::~ContourRecognition()
{

}

QList<QPointF> ContourRecognition::blowUpContour(QList<QPointF> &dataList, double blow, QImage &image)
{
//    QVector2D centroid(0,0);
//    double num = 0;

//    foreach(QPointF current, dataList){
//        centroid.setX(centroid.x()+current.x());
//        centroid.setY(centroid.y()+current.y());
//        num++;
//    }
//    centroid.setX(centroid.x() / num);
//    centroid.setY(centroid.y() /num);

    QList<QPointF> blownList;

    QVector3D normal1, normal2;
    for(int i = 0; i < dataList.size()-1; i++){
        QVector3D halfPoint((dataList.at(i).x()+dataList.at(i+1).x())/2.0, (dataList.at(i).y()+dataList.at(i+1).y())/2.0, 0.0);
        normal1.setX(dataList.at(i+1).y() - dataList.at(i).y());
        normal1.setY(-1*(dataList.at(i+1).x() - dataList.at(i).x()));
        normal1.setZ(0);
        normal1.normalize();
        normal2.setX(-1*(dataList.at(i+1).y() - dataList.at(i).y()));
        normal2.setY(dataList.at(i+1).x() - dataList.at(i).x());
        normal2.setZ(0);
        normal2.normalize();

        QVector3D distVector(dataList.at(i+1).x() - dataList.at(i).x(), dataList.at(i+1).y() - dataList.at(i).y(), 0);
        QVector3D scalar1 = QVector3D::crossProduct(distVector, normal1).normalized();
        QVector3D scalar2 = QVector3D::crossProduct(distVector, normal2).normalized();

        if(scalar1.z() < 0)
            blownList.append(QPointF(halfPoint.x() + blow*normal1.x(), halfPoint.y() + blow*normal1.y()));
        if(scalar2.z() < 0)
            blownList.append(QPointF(halfPoint.x() + blow*normal2.x(), halfPoint.y() + blow*normal2.y()));
    }

//    QVector3D intersect1, intersect2;

//    for(int i = 0; i < dataList.size()-1; i++){
//        double dist = (QVector2D(dataList.at(i+1).x(), dataList.at(i+1).y()) - QVector2D(dataList.at(i).x(), dataList.at(i).y())).length();
//        if(sqrt(4*blow*blow/(dist*dist) - 1) > 0){
//            intersect1.setX(1/2.0*(dataList.at(i).x()+dataList.at(i+1).x()) + 1/2.0*sqrt((4*blow*blow)/(dist*dist) - 1)*(dataList.at(i+1).y()-dataList.at(i).y()));
//            intersect1.setY(1/2.0*(dataList.at(i).y()+dataList.at(i+1).y()) + 1/2.0*sqrt((4*blow*blow)/(dist*dist) - 1)*(dataList.at(i).x()-dataList.at(i+1).x()));
//            intersect1.setZ(0);
//            intersect2.setX(1/2.0*(dataList.at(i).x()+dataList.at(i+1).x()) - 1/2.0*sqrt((4*blow*blow)/(dist*dist) - 1)*(dataList.at(i+1).y()-dataList.at(i).y()));
//            intersect2.setY(1/2.0*(dataList.at(i).y()+dataList.at(i+1).y()) - 1/2.0*sqrt((4*blow*blow)/(dist*dist) - 1)*(dataList.at(i).x()-dataList.at(i+1).x()));
//            intersect2.setZ(0);

//            QVector3D distVector(dataList.at(i+1).x() - dataList.at(i).x(), dataList.at(i+1).y() - dataList.at(i).y(), 0);

//            QVector3D scalar1 = QVector3D::crossProduct(distVector, intersect1).normalized();
//            QVector3D scalar2 = QVector3D::crossProduct(distVector, intersect2).normalized();

//            if((centroid - intersect1).length() >= (centroid - intersect2).length())
//                blownList.append(QPointF(intersect1.x(), intersect1.y()));
//            else
//                blownList.append(QPointF(intersect2.x(), intersect2.y()));
//        }
//    }

    return blownList;
}

void ContourRecognition::skeletonize(int posX, int posY, QImage &image)
{
    double num = 0;
    double xPos = 0, yPos = 0;

    for(int i = posX-1; i <= posX+1; i++){
        for(int j = posY-1; j <= posY+1; j++){
            if(image.pixelColor(i,j) != Qt::white){
                xPos += i;
                yPos += j;
                num++;
            }
            image.setPixelColor(i,j, Qt::white);
        }
    }

    if(num != 0){
        xPos /= num;
        yPos /= num;
        QPoint result(qRound(xPos), qRound(yPos));
        image.setPixelColor(result, Qt::black);
    }
}

QList<QPointF> ContourRecognition::scaleDataToImage(QString dataPath, QImage &image)
{
    QList<QPointF> dataList;
    QList<double> xList, yList;
    QFile file(dataPath);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug("Opening error.");
    }
    double x, y, xmin, xmax, ymin, ymax;
    QTextStream stream(&file);
    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        linestream >> x >> y;
        dataList.append(QPointF(x, y));
        xList.append(x);
        yList.append(y);
    }

    file.close();

    qSort(xList);
    qSort(yList);
    xmin = xList.first();
    xmax = xList.last();
    ymin = yList.first();
    ymax = yList.last();

    QList<QPointF> scaledList;

    foreach(QPointF current, dataList){
        QPointF scaled;
        scaled.setX((image.width()/qAbs(xmax-xmin))*(current.x() - xmin));
        scaled.setY(-1*(image.height()/qAbs(ymax-ymin))*(current.y() - ymax));
        scaledList.append(scaled);
    }

//    qDebug("%f %f %f %f", xmax-xmin, ymax-ymin, xmin, ymax);
    return scaledList;
}

QList<QPointF> ContourRecognition::scaleContour(QString dataPath, QImage &image)
{
    QList<QPointF> dataList;
    QList<double> xList, yList;
    QFile file(dataPath);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug("Opening error.");
    }
    double x, y;
    QTextStream stream(&file);
    while(!stream.atEnd()){
        QString line = stream.readLine();
        QTextStream linestream(&line);
        linestream >> x >> y;
        dataList.append(QPointF((image.width()/86.98)*(x - (-67.989)), -1*(image.height()/33.59)*(y - 83.599)));
    }

    file.close();
    return dataList;
}

double ContourRecognition::blowDistance(double R, double s, double H, double h)
{
    double distance;
    double phi, theta;

    phi = acos(R/(R+H));
    theta = acos(R/(R+h));

    distance = R*(phi+theta) - s;

    return distance;
}

QList<QPointF> ContourRecognition::getRelativeContourPositions(QImage &image)
{
    QList<QPointF> relContPos;
    for(int i = 0; i < image.width(); i++){
        for(int j = 0; j < image.height(); j++){
            if(image.pixelColor(i,j) == Qt::blue)
                relContPos.append(QPointF((double)i/(double)image.width(), (double)j/(double)image.height()));
        }
    }
    return relContPos;
}

