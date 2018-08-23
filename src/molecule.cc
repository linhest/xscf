#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <vector>

#include "molecule.h"
#include "atom_labels.h"
using namespace std;

extern const char atom_labels[][4];

#define MAXLINELENGTH 2000
molecule::molecule (char filename[]) {

  // default values
  energy_conv = 1e-8;
  density_conv = 1e-8;
  n_diis = 8;

  print = 0;
  filename_guess_save[0] = '\0';
  filename_guess_read[0] = '\0';

  ifstream fin;
  fin.open (filename);
  if (!fin.good ()) {
    throw std::invalid_argument (strcat (filename, " not found "));
  }
  while (!fin.eof ()) {
    char buf[MAXLINELENGTH];
    fin.getline (buf, MAXLINELENGTH);

    int l = strlen (buf);
    if (l == MAXLINELENGTH) {
      throw std::invalid_argument ("line too long");
    }
    if (strchr (buf, '#') != NULL) {
      char * a = strchr (buf, '#');
      *a = '\0';
    }
    // concatenate lines if end with \\n
    while (l > 2 && buf[l - 1] == '\\') {
      fin.getline (buf + (l - 1), MAXLINELENGTH - l);
      l = strlen (buf);
      if (l == MAXLINELENGTH)
        throw std::invalid_argument ("line too long");
    }

    char * token = strtok (buf, "=");
    if (token == NULL)
      continue;
    char * token2 = strtok (NULL, "=");
    if (token2 == NULL)
      continue;
    if (strcasecmp (token, "nref") == 0) {
      nref = parse_int (token2);
    } else if (strcasecmp (token, "print") == 0) {
      print = parse_int (token2);
    } else if (strcasecmp (token, "basis_set") == 0) {
      strncpy (basis_set, token2, LEN_FILENAME);
    } else if (strcasecmp (token, "guess_save") == 0) {
      strncpy (filename_guess_save, token2, LEN_FILENAME);
    } else if (strcasecmp (token, "guess_read") == 0) {
      strncpy (filename_guess_read, token2, LEN_FILENAME);
    } else if (strcasecmp (token, "exlevels") == 0) {
      exlevels = parse_int (token2);
    } else if (strcasecmp (token, "occ") == 0) {
      locc=count_entries(token2, ',');
      occ = new unsigned int[locc];
      parse_int_array_fixed_len (occ, locc, token2, ',');
    } else if (strcasecmp (token, "charge") == 0) {
      charge = parse_int (token2);
    } else if (strcasecmp (token, "sym") == 0) {
      target_sym = parse_int (token2);
    }  else if (strcasecmp (token, "nmo") == 0) {
      //      nmo = parse_int (token2);
    } else if (strcasecmp (token, "ref") == 0) {
      // vector<char *>* list = parse_string_array (token2, ';');
      // nref = list->size ();
      // int n=count_entries((*list)[0], ',');
      // refs = new unsigned int[nref * n];
      // int i;
      // for (i = 0; i < nref; i++) {
      //   parse_int_array_fixed_len (&refs[i * n], n, (*list)[i], ',');
      // }
      // delete list;
    } else if (strcasecmp (token, "geom") == 0) {
      char * ptr_start = token2;
      l = strlen (token) + strlen (token2);
      ptr_start = strchr (token2, '(') - 1;
      while (ptr_start == NULL) {
        fin.getline (buf + l - 1, MAXLINELENGTH);
        ptr_start = strchr (ptr_start, '(');
        l = strlen (token) + strlen (token2);
      }
      // concatenate lines until ')'
      while (strchr (ptr_start, ')') == NULL) {
        *(ptr_start + strlen (ptr_start)) = '\n';
        fin.getline (ptr_start + strlen (ptr_start) + 1,
        MAXLINELENGTH - strlen (token) - strlen (token2));
        l = strlen (token) + strlen (token2) + 3;
        if (l >= MAXLINELENGTH)
          throw std::invalid_argument ("line too long");
      }
      char * ptr_end = strchr (ptr_start, ')');
      *(ptr_end + 1) = '\0';
      *(ptr_end) = '\0';
      ptr_start = ptr_start + 2;
      l = strlen (buf);
      if (l == MAXLINELENGTH)
        throw std::invalid_argument ("line too long");
      vector<double*> G_list;
      vector<char*> label_list;
      int i;
      char *ptr, *ptr_next;
      for (i = 0, ptr = ptr_start + 1; ptr != NULL;) {
        char label[4];
        double g[3];
        int num;
        ptr_next = strchr (ptr, '\n');
        if (ptr_next != NULL)
          *(ptr_next) = '\0';
        if (strlen (ptr) != 0) {
          num = sscanf (ptr, "%s %lf %lf %lf ", label, &g[0], &g[1], &g[2]);
          if (num == 4) {
            double * g_new = new double[3];
            char * l_new = new char[4];
            g_new[0] = g[0];
            g_new[1] = g[1];
            g_new[2] = g[2];
            strncpy (l_new, label, 4);
            G_list.push_back (g_new);
            label_list.push_back (l_new);
            i++;
          } else if (num > 0) {
            throw std::invalid_argument ("error reading geom");
          }
        }
        if (ptr_next == NULL)
          break;
        ptr = ptr_next + 1;
      }
      natm = G_list.size ();
      geom = new double[natm * 3];
      Z = new int[natm];
      atm_label = new char *[natm];
      for (i = 0; i < natm; i++) {
        geom[i * 3 + 0] = G_list[i][0];
        geom[i * 3 + 1] = G_list[i][1];
        geom[i * 3 + 2] = G_list[i][2];
        atm_label[i] = label_list[i];
        delete[] G_list[i];

        for (int j = 0; j < N_ATOM_LABELS; j++) {
          if (strncasecmp (atm_label[i], atom_labels[j], 4) == 0)
            Z[i] = j;
        }
      }
      G_list.empty ();
      label_list.empty ();
    } else if (strcasecmp (token, "energy_conv") == 0) {
      energy_conv = parse_float (token2);
    } else if (strcasecmp (token, "density_conv") == 0) {
      density_conv = parse_float (token2);
    } else {
      throw std::invalid_argument (strcat (token, " not recognized"));
    } // check token
  } //!feof
  fin.close ();

  Enuc = calculate_Enuc (geom, Z, natm);

} // initialize

molecule::~molecule () {
  for (int i = 0; i < natm; i++) {
    delete[] atm_label[i];
  }
  delete[] atm_label;
  //  delete[] refs;
  delete[] Z;
  delete[] geom;
  delete[] occ;

}

ostream&
operator<< (ostream &out, const molecule & mol) {
  int i;
  //  out << "#nmo=" << mol.nmo << endl;
  out << "#ms2=" << mol.ms2 << endl;
  out << "#nref=" << mol.nref << endl;
  out << "#basis_set=" << mol.basis_set << endl;
  out << "#Geometry=" << endl;
  for (i = 0; i < mol.natm; i++) {
    printf ("# %d %d %4s %f %f %f\n", i, mol.Z[i], mol.atm_label[i],
            mol.geom[i * 3 + 0], mol.geom[i * 3 + 1], mol.geom[i * 3 + 2]);
  }
  out << "#Enuc=" << mol.Enuc << endl;

  // out << "#References=" << endl;
  // for (i = 0; i < mol.nref; i++) {
  //   out << "#" << i << ":";
  //   for (j = 0; j < mol.nmo; j++) {
  //     out << mol.refs[i * mol.nmo + j];
  //     if (j < mol.nmo - 1) {
  //       out << ",";
  //     }
  //   }
  //   out << endl;
  // }
  return (out);
}

int
molecule::get_natm () {
  return (natm);
}

int
molecule::parse_int (char * string) {
  int n;
  stringstream ss;
  ss << string;
  ss >> n;
  return (n);
}

double
molecule::parse_float (char * string) {
  double f;
  stringstream ss;
  ss << string;
  ss >> f;
  return (f);
}

vector<char *>*
molecule::parse_string_array (char * buf, const char delim = ',') {
  vector<char *>* s = new vector<char *> ();
  char * token;
  for (token = strtok (buf, &delim); token != NULL;
      token = strtok (NULL, &delim)) {
    s->push_back (token);
  }
  return (s);
}

unsigned int *
molecule::parse_int_array_fixed_len (unsigned int * s, int n, char * buf,
                                     const char delim = ',') {
  char * token;
  int i;
  for (i = 0, token = strtok (buf, &delim); token != NULL && i < n; token =
      strtok (NULL, &delim), i++) {
    s[i] = parse_int (token);
  }
  if (i == n && token != NULL)
    throw std::invalid_argument ("array too long");
  for (; i < n; i++) {
    s[i] = 0;
  }

  return (s);
}

int molecule::count_entries (char * buf, const char delim = ',') {
  int i;
  char * ptr;
  for (i = 1, ptr = strchr (buf, (int)delim); ptr != NULL; i++){
    ptr=strchr(ptr+1,(int)delim);
  }
  return(i);
}

int *
molecule::parse_int_array_fixed_len (int * s, int n, char * buf,
                                     const char delim = ',') {
  char * token;
  int i;
  for (i = 0, token = strtok (buf, &delim); token != NULL && i < n; token =
      strtok (NULL, &delim), i++) {
    s[i] = parse_int (token);
  }
  if (i == n && token != NULL)
    throw std::invalid_argument ("array too long");
  for (; i < n; i++) {
    s[i] = 0;
  }
  return (s);
}

double
calculate_Enuc (double * geom, int * Z, int natm) {
  double E = 0.;
  for (int i = 0; i < natm; i++) {
    for (int j = i + 1; j < natm; j++) {
      double dx = (geom[i * 3 + 0] - geom[j * 3 + 0]);
      double dy = (geom[i * 3 + 1] - geom[j * 3 + 1]);
      double dz = (geom[i * 3 + 2] - geom[j * 3 + 2]);
      double dist = sqrt (dx * dx + dy * dy + dz * dz);
      cout << i << " " << j << " " << dist << endl;
      E = E + Z[i] * Z[j] / dist;
    }
  }
  return (E);
}

