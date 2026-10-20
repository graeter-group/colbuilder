# Glucosepane: ottimizzazione sterica e diagnostica dei contatti

## Stato del fix

Implementazione locale del 2026-09-06, attiva per i crosslink AGS/LGX nella
sequence generation. Non modifica le parametrizzazioni, i termini bonded,
le pairs 1-4 o le exclusions. Le mosse di ottimizzazione PYD e MOLD restano
quelle preesistenti. Non sono necessarie nuove chiavi YAML.

La prima revisione correggeva il clash specifico di 30_bis_GCP sul sito
310.C - 542.B, ma bloccava anche 523.A - 286.C: i suoi limiti sterici erano
troppo restrittivi per un filtro pre-minimizzazione. La revisione successiva
separava le penalita' di ottimizzazione e gli avvisi dai soli contatti quasi
coincidenti che impedivano l'esportazione. Su richiesta dell'utente, il controllo
finale del GRO e' ora solo diagnostico: tutti i contatti con coordinate finite,
anche estremi o esattamente coincidenti, producono warning senza bloccare
l'esportazione. Restano errori le coordinate non finite e le incoerenze GRO/ITP.
La ricerca in sequence generation mantiene invece i propri criteri di
accettazione; questa modifica riguarda soltanto il controllo finale sul GRO.
Le soglie condivise sono in `core/utils/steric_policy.py`; non sono parametri
del force field ne' un certificato di stabilita' MD.

## Prima

La vecchia MC accettava i candidati in base alla distanza forming, senza un
controllo sterico dell'ambiente. Alcune mosse ruotavano arbitrariamente la
catena laterale rispetto al backbone o intervenivano su porzioni del backbone.
Il successo sulla distanza poteva quindi accompagnarsi a sovrapposizioni
catastrofiche con altre catene della stessa tripla elica.

Nel caso originale, AGS310.C:NE si avvicinava a GLU312.A:CB da 6.54 A a 0.93 A.
Dopo pdb2gmx, NE e HB1 erano a circa 0.124 A: la minimizzazione falliva con
forza infinita e LJ(SR) dell'ordine di 2.9e17 kJ/mol.

## Nuova ricerca geometrica

File: `src/colbuilder/core/sequence/steric_optimization.py`.

1. Il grafo intraresiduo e' letto dal file RTP AMBER distribuito con Colbuilder.
   Si aggiungono i legami peptidici e il legame futuro AGS:NZ - LGX:CE.
   Le esclusioni del punteggio sono basate sul grafo, non sulle distanze.
2. Le mosse sono torsioni attorno a legami aciclici della catena laterale che
   coinvolgono un carbonio alifatico CT. Non si tagliano anelli o legami
   coniugati; non si muovono N, CA, C, O o gli idrogeni del backbone.
3. Queste rotazioni conservano matematicamente lunghezze e angoli dei legami
   interni e muovono rigidamente gli anelli. La lettura/scrittura PDB e le copie
   Chimera introducono il normale arrotondamento delle coordinate.
4. Le occorrenze dello stesso residuo nelle due copie ricevono la stessa
   trasformazione. Questo e' coerente con il seed che verra' replicato.
5. Il punteggio comprende chiusura del legame, sovrapposizione sterica soft,
   una penalita' geometrica per gli angoli del nuovo legame e penalita' forti
   per la chiusura desiderata e i contatti quasi coincidenti. I contatti
   compressi non ricevono piu' la barriera forte del vecchio limite di 2 A.
6. La ricerca usa `scipy.optimize.differential_evolution`, seguita dal polishing
   locale di SciPy. Il budget per ricerca e' 500 generazioni, popolazione di
   10 volte il numero di torsioni. Il seed dipende dallo stato NumPy.
7. Se il peggior ostacolo e' una catena laterale standard con torsioni disponibili,
   possono essere aggiunti fino a due residui vicini alla ricerca. Anche per
   loro backbone, legami e angoli interni restano invariati. Non si modificano
   altri marker di crosslink con questa espansione.
8. Si conserva il miglior candidato ammissibile incontrato, anche se l'ultimo
   punto restituito dall'ottimizzatore non fosse ammissibile.
9. Il candidato viene ricontrollato dopo il riporto delle coordinate sulle
   copie e sul seed. In caso di errore, il riporto viene annullato.

Il punteggio sterico usa raggi elementari e penalita' quadratiche. Le coppie
1-2 e 1-3 non ricevono il punteggio non legato; le 1-4 sono incluse con un
trattamento geometrico piu' morbido. Questo NON e' il calcolo AMBER delle
energie: non sostituisce LJ, elettrostatica, termini torsionali o minimizzazione.
Le topologie conservano invece le loro interazioni e parametrizzazioni esatte.

## Criteri di accettazione della sequenza

- Forming AGS:NZ - LGX:CE: obiettivo di ricerca 1.2-2.5 A; 2.5-5.0 A produce
  un avviso, non un errore sterico. Restano errori i forming sotto 1.2 A o oltre
  5.0 A. Il limite massimo deriva da `MAX_DIVALENT_DISTANCE`, gia' pari a 5 A
  in Colbuilder e coerente con il raggruppamento atomistico dei crosslink.
  La vecchia annotazione che chiamava 2.5 A il limite preesistente era errata.
- Contatto non legato heavy-heavy sotto 0.25 A: errore, anche per le coppie 1-4.
- Contatto heavy-heavy ordinario fra 0.25 e 2.0 A: accettato con avviso.
- Contatto heavy-heavy 1-4 fra 0.25 e 1.5 A: accettato con avviso.
- Punteggio e coordinate devono restare validi e finiti. Le altre verifiche
  strutturali (atomi mancanti, copie incoerenti, backbone) non sono disabilitate.

Il punteggio continua a favorire legami vicini a 1.5 A e a penalizzare tutti i
contatti compressi: rendere un contatto non bloccante non lo rende desiderabile.
Si cerca un margine di 0.05 A rispetto alla soglia catastrofica. Una struttura
accettata e' soltanto un candidato pre-minimizzazione, non una struttura
rilassata o certificata MD. 0.25/0.15 A sono soglie pragmatiche di quasi
coincidenza, NON una previsione universale del limite delle forze in GROMACS.

Il caso 10_GCP salvato dall'utente contiene AGS:C16 - LGX:CE a circa 0.323 A.
Il cammino e' C16-C15-NZ-CE: e' una coppia 1-4 esplicitamente presente nella
topologia. L'utente riporta minimizzazione e NVT riuscite per questo sistema;
la revisione non equipara quindi tale contatto a una forza necessariamente
infinita. Lo segnala comunque come compresso. Non viene tolta alcuna pair.

## Integrazione e diagnostica

- `sequence/optimize_crosslinks.py`: solo AGS/LGX passa al nuovo percorso;
  le modifiche sono mantenute nel seed e tutte le ottimizzazioni GCP della
  chiamata vengono ricontrollate prima dell'esportazione.
- `utils/crosslinks.py`: un candidato non ammissibile genera
  `StericOptimizationError`. Il wrapper riprova entro il budget gia' esistente
  di tre tentativi e infine solleva `SequenceGenerationError` senza esportare
  un nuovo PDB dichiarato valido. Gli avvisi di compressione non attivano
  questo percorso di errore e non consumano ulteriori tentativi.
- `sequence/sequence_generator.py`: gli errori scientifici gia' classificati
  vengono propagati invece di trasformarli in un errore generico di setup
  delle directory. I log riportano distanza forming e contatto peggiore.

## Controllo sul GRO completo

File: `topology/coordinate_validation.py`, chiamato da `build_amber99()` dopo
la scrittura del GRO, prima di completare con successo l'esportazione.

Il controllo legge gli ITP nell'ordine del GRO, verifica identita' e conteggio
degli atomi e ricostruisce i legami. Cerca contatti gravi che coinvolgono AGS
o LGX, includendo idrogeni e contatti fra molecole/ITP diversi:

- heavy-heavy sotto 0.25 A: warning severo, senza blocco;
- contatti con almeno un idrogeno sotto 0.15 A: warning severo, senza blocco;
- heavy-heavy sotto 2.0 A (1.5 A per le 1-4), oppure con H sotto 0.6 A,
  ma sopra le soglie severe: avviso, con esportazione consentita;
- legami e coppie 1-3 non vengono trattati come contatti non legati.

Il log elenca TUTTE le coppie sotto le soglie severe, ordinate per distanza,
con residuo e numero residuo, nome atomo, tipo del force field, indice globale
1-based nell'ordine del GRO, nome ITP e indice locale 1-based. L'indice globale
non e' il seriale GRO che ricomincia da zero dopo 99999. Il numero residuo
proviene dall'ITP: non e' un identificativo di catena; ITP e indice locale
identificano senza ambiguita' l'atomo anche con numeri residuo ripetuti.
Per i contatti meno estremi vengono riportati conteggio e coppia peggiore.

La funzione restituisce i conteggi reali in `gross_overlaps` e
`compressed_contacts`; il successo della chiamata non significa assenza di
clash. I due conteggi sono disgiunti. I messaggi compaiono in console e nel
log Colbuilder tramite il logger standard. Il controllo non sposta atomi,
non cambia parametri, non rimuove interazioni e non aggiunge exclusions.

Il vecchio GRO di 30_bis_GCP viene ora esportato con 60 warning severi,
compreso il contatto HB1-NE di circa 0.124 A: resta un caso noto di forza
infinita, NON e' stato corretto dal cambio di policy. Il GRO salvato di 10_GCP
continua a essere accettato con avvisi.

La motivazione e' il test della terna originale rigenerata: dieci contatti
AGS494:HO18 - GLY266:C a circa 0.123 A bloccavano il filtro geometrico, ma la
minimizzazione diagnostica ha completato 500 passi con forze finite. Il tipo
HO non ha Lennard-Jones in questo force field. Dopo EM, le distanze erano
3.66-3.78 A e non restavano warning GCP. Non si puo' dedurre la stessa
stabilita' per qualunque contatto sotto soglia: controllare sempre energia,
forze finite e convergenza nella preparazione della simulazione.

## Limiti importanti

La ricerca torsionale vede le due triple eliche complete richieste da
crosslink_copies, non tutte le immagini del cristallo. Il controllo sul GRO
successivo amplia la verifica a tutta la fibrilla effettivamente esportata,
ma usa coordinate nonperiodiche ed e' una diagnostica non bloccante.
Le condizioni periodiche della box MD finale, gli idrogeni aggiunti in fasi
successive, il solvente e la convergenza energetica vanno validati separatamente.

Il controllo globale riguarda GCP; non certifica la geometria PYD/MOLD o dei
capping groups. I vecchi clash gia' presenti in questi gruppi non vengono
corretti automaticamente. In particolare, una distanza oltre la soglia minima
di questo filtro non equivale a un contatto privo di tensione sterica.

La revisione delle soglie o del comportamento di export non costituisce una
nuova validazione MD/NVT/pulling.
La terna alternativa studiata in `overlap_search_20260906` aveva superato 500
passi di minimizzazione diagnostica con la policy precedente, senza raggiungere
la convergenza. Quel risultato non certifica nuovi PDB prodotti da questa
revisione e non sostituisce la preparazione del sistema simulato.

## Test e riproduzione

- `tests/test_glucosepane_steric_optimization.py`: invarianza di backbone,
  legami, angoli e anelli; maschere dal grafo; copie coerenti; rifiuto di
  candidati patologici; rollback; conservazione del miglior punto ammissibile.
- `tests/test_glucosepane_coordinate_validation.py`: warning per clash heavy
  e con H, anche coincidenti e fra ITP diversi; metadati dei responsabili,
  conteggi senza duplicazioni, invariabilita' degli input, maschere 1-2/1-3
  e warning 1-4; errori mantenuti per mapping e coordinate non finite.
- `tests/run_saved_glucosepane_optimization.py`: prova su copie salvate prima
  della vecchia MC, senza Chimera; salva PDB e metriche in una directory nuova.
- `tests/run_glucosepane_steric_regression.py`: esegue i due YAML del caso in
  directory isolate, con eventuale geometria/topologia solo dopo il successo
  di entrambi i passaggi di sequenza.
- `tests/run_glucosepane_policy_regression.py`: verifica che entrambi i GRO
  10_GCP e 30_bis_GCP siano esportabili, ma che il secondo mantenga il
  conteggio non nullo dei contatti severi; poi riesegue il YAML 523.A - 286.C, D2-D1,
  in una directory nuova. Verifica successo, contatti e backbone invariato.

Alla verifica precedente della policy del 2026-09-06: 242 test automatici superati.
La dipendenza pytest
di prova e' stata installata in `/private/tmp/colbuilder-test-dependencies`,
senza cambiare le dipendenze dell'ambiente Colbuilder.

Le prove reali sono in:
`/Users/guidogiannetti/Desktop/kimmdy_age_campaign/glucosepane_rattus/30_bis_GCP/steric_fix_test_20260906`.

Per 310.C - 542.B:

| Seed | Forming (A) | Min. non legato heavy (A) | Spostamento massimo backbone (A) |
|---:|---:|---:|---:|
| 1201 | 2.4520 | 2.1193 | 0 |
| 1202 | 2.4528 | 2.0966 | 0 |

Per il seed 1201, NE310.C - CB312.A e' ora 4.950 A. Nessun atomo al di fuori
dei due residui GCP cambia posizione nel PDB; lo scarto massimo nelle lunghezze
dei legami interni dei marker, inclusi gli arrotondamenti PDB/Chimera, e' 0.0043 A.

Questi erano risultati della prima policy: il caso completo si arrestava su
523.A - 286.C anche permettendo torsioni di GLU517.B; un test separato su
494.B - 268.A restituiva forming circa 5.15 A e contatto minimo circa 1.32 A.
Non erano prove di impossibilita' fisica. Il secondo caso resta oltre il
limite di connettivita' di 5 A: una soglia sterica piu' permissiva non garantisce
che ogni sito riesca a chiudersi correttamente.

Risultati della revisione corrente:
`/Users/guidogiannetti/Desktop/kimmdy_age_campaign/glucosepane_rattus/10_GCP/steric_policy_final_20260906`.
Il controllo salva `result.json`, il YAML realmente eseguito e tutti gli
intermedi nella sottocartella `sequence`. I file originali non sono sovrascritti.

Esito storico della regressione prima del passaggio del GRO a warning-only:

- GRO originale 10_GCP: accettato, 300 contatti compressi segnalati come avvisi,
  nessuna quasi-coincidenza oltre le nuove soglie di errore.
- GRO originale 30_bis_GCP: ancora rifiutato, 60 quasi-coincidenze rilevate.
- Nuova sequenza 523.A - 286.C con D2-D1: successo al primo tentativo; forming
  2.4617 A, minimo heavy nonbonded 1.2001 A (avviso), minimo heavy 1-4 2.7222 A,
  spostamento massimo del backbone 0.000 A.
- Non sono state rieseguite minimizzazione o NVT di questo nuovo seed. I test
  sui GRO sopra usano gli output originali salvati, non una nuova topologia.

## Verifica del comportamento warning-only

Revisione finale del 2026-09-06: 253 test automatici superati, incluse 36 prove
del controllo coordinate. I warning di deprecazione preesistenti Pydantic/SWIG
non riguardano questa modifica. `git diff --check` senza errori.

La topologia della terna originale e' stata rigenerata dai caps conservati,
senza rifare la sequenza o cambiare coordinate, mediante
`tests/run_saved_system_regression.py --temporary-topology`. Risultato:

- Esportazione completata, 30 GCP + 20 PYD, 35 gruppi molecolari.
- Dieci coppie severe elencate nel log, piu' 272 contatti compressi nel riepilogo.
- `gmx grompp` completato con codice 0, senza `-maxwarn`.
- Tutti i 104 file TOP/GRO/ITP/force field confrontati sono identici byte per
  byte ai temporanei precedenti: cambia il controllo di flusso, non il sistema.
- I GRO originali 10_GCP e 30_bis_GCP non sollevano piu' errori sterici:
  mantengono rispettivamente 0 e 60 contatti severi nelle metriche restituite.

Output e report riproducibili:
`/Users/guidogiannetti/Desktop/kimmdy_age_campaign/glucosepane_rattus/30_bis_GCP/original_sites_regenerated_20260906/warning_only_export_20260906`.

L'ambiente conda `colbuilder` importa gia' questo checkout in modalita'
editable; non occorre reinstallarlo ne' cambiare i file YAML. Nessuna nuova
MD/NVT/pulling e' stata eseguita durante questa modifica.
