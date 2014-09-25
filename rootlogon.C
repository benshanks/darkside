{
char *inc = gSystem->ExpandPathName("$CLHEP_INCLUDE_DIR");
gInterpreter->AddIncludePath(inc);
delete [] inc;
}