# Ratio replacement: modifiche al codice, passo per passo

Questa nota descrive la correzione della continuita del backbone introdotta il
6 settembre 2026. Riguarda la mutazione dei marker in residui standard, non la
scelta di quali crosslink rimuovere. Il documento tecnico in inglese e
[backbone_preservation.md](backbone_preservation.md).

Aggiornamento: la correzione successiva del conteggio del ratio e del grafo dei
crosslink finali e descritta in [crosslink_network.md](crosslink_network.md).
La presente nota riguarda il fix del backbone, non tutte le modifiche successive.

## 1. Il difetto originale

Nel precedente `pyd_50pp` cinque legami peptidici LYS103:C--GLY104:N erano
assenti dagli ITP. Il PDB prima del replacement conteneva un C-N allungato,
circa 2.30 A, ma nessuna terminazione tra quei residui.

Durante `open` / `swapaa` / `write`, Chimera poteva reinterpretare quel contatto
come una separazione di catena e introdurre un `TER`. Pdb2gmx leggeva quindi
due estremita e produceva CLYS/NGLY, anziche il legame peptidico interno.
Il vecchio `repair_missing_backbone_atoms` ripristinava atomi mancanti, ma non
ricostruiva i confini polimerici originali: da solo non risolveva il problema.

Non bastava aggiungere il bond al vecchio ITP: terminazioni artificiali
comportano anche atomi, cariche, tipi e termini bonded differenti.

## 2. Nuovo modulo e lettura dell'input

File: [core/geometry/backbone.py](../src/colbuilder/core/geometry/backbone.py).

`read_pdb_residues` legge i record PDB mantenendo l'ordine e senza perdere
`TER`, cambi di chain ID e confini MODEL/ENDMDL. Ogni residuo conserva:

- identita: modello, chain ID, numero residuo e insertion code;
- nome del residuo;
- identificatore del segmento polimerico;
- record originali degli atomi e loro posizione nel file.

La continuita e definita dai segmenti dell'input, non da una soglia spaziale.
Due residui entro lo stesso segmento restano consecutivi anche se il loro
C-N e troppo lungo. Una vera terminazione resta tale anche se le due parti
sono vicine o usano lo stesso chain ID. Nessuna numerazione mancante viene
usata per inventare residui intermedi.

Il parser rifiuta record troncati, atomi duplicati e alternate locations
ambigue. Prima della mutazione si rifiutano anche piu modelli o identita di
residuo ripetute, perche il comando swapaa non le indirizzerebbe univocamente.

## 3. Validazione delle istruzioni e snapshot

`preserve_backbone(type_dir, instructions)` e un context manager.
Prima di chiamare Chimera:

1. Legge istruzioni `file.caps.pdb LYS/ARG numero chain`.
2. Verifica che il file rimanga nella directory consentita.
3. Risolve un unico residuo, con matching esatto della chain quando possibile.
4. Rifiuta target mancanti, istruzioni incompatibili e target diversi da LYS/ARG.
5. Richiede N, CA, C e O nel backbone originale del residuo da mutare.
6. Salva in memoria il contenuto completo di tutti i file coinvolti.

Lo snapshot avviene prima di qualsiasi mutazione del batch. Non vengono
modificati ne `ratio_replace` ne la selezione casuale ne il numero di istruzioni.

## 4. Chimera resta responsabile delle nuove catene laterali

Non e stato modificato `chimera_scripts/swapaa.py` da questo intervento.
L'invocazione esistente viene racchiusa nella protezione:

```python
with preserve_backbone(type_dir, instructions):
    success = await self._run_chimera_command(...)
    if not success:
        raise GeometryGenerationError(...)
```

Il PDB completo riscritto da Chimera non viene accettato come nuova autorita
sulla connettivita: da esso si ricavano soltanto le catene laterali richieste.

## 5. Ricostruzione del PDB dopo swapaa

`_restore_sidechains` confronta l'output Chimera con lo snapshot:

- richiede gli stessi residui nello stesso ordine e con la stessa identita;
- richiede il nome LYS/ARG atteso nei target e nessuna mutazione altrove;
- richiede gli atomi pesanti completi della nuova catena laterale;
- controlla che le coordinate degli atomi laterali siano finite;
- rifiuta spostamenti maggiori di 0.01 A degli N/CA/C/O presenti nell'output
  del target, per non innestare una catena laterale da un riferimento spostato.

Costruisce poi il risultato usando il PDB originale come struttura di base:

1. I residui non mutati provengono interamente dall'originale.
2. N/CA/C/O e gli altri atomi di backbone/terminazione riconosciuti dei target
   provengono dall'originale, cambiando soltanto il nome del residuo.
3. Le vecchie catene laterali dei marker sono rimosse e sostituite con quelle
   costruite da Chimera.
4. Sono mantenuti solo i confini polimerici originali. I nuovi TER di Chimera
   non entrano nel risultato; i TER autentici e i caps originali rimangono.

Un ossigeno O perso da Chimera viene quindi recuperato dallo snapshot.
Le coordinate finali del backbone non sono approssimate o ottimizzate: sono
quelle originali, anche sotto la tolleranza di 0.01 A.

I seriali ATOM/HETATM/TER vengono rinumerati. I CONECT relativi ad atomi
conservati vengono rimappati; gli archi che coinvolgevano catene laterali
sostituite sono eliminati. La chimica delle nuove catene laterali viene dai
template RTP, non dai vecchi CONECT. ANISOU e MASTER sono gestiti evitando
riferimenti o conteggi divenuti obsoleti.

## 6. Errore e ripristino dei file

Tutti i risultati vengono validati prima di iniziare a pubblicarli.
`_atomic_write` scrive un file temporaneo e usa `os.replace` per sostituire
ciascun PDB. Se Chimera, la validazione o la scrittura falliscono, la protezione
tenta di ripristinare tutti i file interessati dai byte salvati nello snapshot,
anche quando Chimera ha cancellato un file.

La gestione include le cancellazioni Python. Non e una transazione garantita
contro spegnimenti, kill del processo o guasti del disco: l'atomicita vale per
il singolo file e un guasto che impedisce le scritture puo impedire il rollback.

## 7. Tutti i punti di integrazione

In [geometry_replacer.py](../src/colbuilder/core/geometry/geometry_replacer.py):

- `replace_in_system`: protegge il replacement sul System, incluso il ratio.
- `replace_direct`: protegge il replacement su PDB fornito direttamente.
- `_apply_manual_replacements_to_dir`: protegge istruzioni manuali e auto-fix.

In [crosslink_mixer.py](../src/colbuilder/core/geometry/crosslink_mixer.py):

- `_apply_chimera_swaps`: applica la stessa protezione alle mutazioni del mixer.

I vecchi helper atom-only restano disponibili per compatibilita, ma questi
percorsi interni non li usano piu. La protezione riguarda LYS e ARG e non
dipende dallo scope enzymatic/non_enzymatic.

## 8. Seconda barriera: controllo dell'ITP

In [amber.py](../src/colbuilder/core/topology/amber.py), `Amber.write_itp`
chiama `validate_backbone_topology` dopo l'aggiunta dei termini dei crosslink.
Il validatore confronta il PDB effettivamente passato a pdb2gmx con l'ITP:

1. Controlla numero e ordine dei residui mediante i numeri residuo.
2. Richiede C(i)--N(i+1) per le adiacenze previste nello stesso segmento PDB.
3. Richiede N--CA, CA--C e C--O presenti nell'input.
4. Consente O rinominato OC1/O1 soltanto alle vere estremita del segmento.
5. Rifiuta legami peptidici attraverso confini polimerici originali.

Gli indici numerici dell'ITP vengono risolti nei nomi atomici: non si assume
che gli atomi abbiano gli stessi seriali del PDB.

Una violazione diventa `TopologyGenerationError`, codice `TOP_ERR_005`.
`build_amber99` gia rilanciava questa classe di errore: la nuova validazione
usa quel percorso per interrompere la generazione, non saltare il gruppo e
produrre un risultato apparentemente riuscito. Questo controllo e specifico
del backbone AMBER all-atom; non certifica la topologia dei crosslink o Martini.

## 9. Configurazione, test e limiti

Nessun nuovo campo YAML. `contact_distance=25 A`, parametri dei crosslink,
esclusioni, angoli, dihedrals e mosse MC non sono stati modificati dal fix.
L'installazione conda locale e editable e carica direttamente questo repository.

[test_backbone_preservation.py](../tests/test_backbone_preservation.py) aggiunge
30 casi ai 12 test esistenti: TER veri/artificiali, coordinate, O mancante,
LYS/ARG, caps, seriali/CONECT, input ambigui, rollback e validazione AMBER.
Il replay reale riproducibile e
[run_backbone_replacement_regression.py](../tests/run_backbone_replacement_regression.py).

Per vecchi output difettosi occorre ripartire dai PDB **prima del replacement**.
Un TER artificiale gia presente nell'input e indistinguibile da uno intenzionale
senza ulteriore provenienza; questo fix non lo elimina arbitrariamente.

Il controllo preserva e verifica la continuita prevista, non certifica che il
PDB di partenza sia chimicamente corretto. Non corregge distanze anomale,
clash o distorsioni HYP e non dimostra stabilita in MD/pulling. La presenza
dei PYD previsti, l'assenza di marker orfani e il rapporto di replacement vanno
verificati separatamente: un backbone continuo non garantisce questi aspetti.
