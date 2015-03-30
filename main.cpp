#include "simgast.h"
#include <QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	SimGast w;
	w.show();
	return a.exec();
}