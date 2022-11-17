# pythf
script that have some auxiliary function on the post- and pre- processing of TOUGH+HYDRATE simulation  

At current, it can only be used on the MESH with properties shown as follow: 1):2D MESH, 2):structured grid mesh

It can modify the MESH file to change the subdomain into what you set, a sample is shown as follow.
![ini_Shyd](https://user-images.githubusercontent.com/111716358/202487423-38616b62-9762-4fe7-98e9-9445ea2e1579.jpg)

It can caculate the initial pressure and temperature that have distribution with gredient for the active element and generate a INCON file automatically
[INCON.txt](https://github.com/tjuCJC/pythf/files/10033064/INCON.txt)


It can visualize the initial condition in the INCON file and final condition in the SAVE, for example:
![ini_P](https://user-images.githubusercontent.com/111716358/202487674-e74650b7-45d4-4284-94e1-f5b3fb76e275.jpg)

![ini_T](https://user-images.githubusercontent.com/111716358/202487750-e63cae53-8ecd-4765-8011-12b3ff7e1897.jpg)

and
![save_P](https://user-images.githubusercontent.com/111716358/202489062-cdecd401-8dec-447f-bae7-7c609e79a696.jpg)
![save_T](https://user-images.githubusercontent.com/111716358/202489107-3f671b12-3fd0-4fbb-b336-7c441466a18c.jpg)
![save_Sgas](https://user-images.githubusercontent.com/111716358/202489143-ca45c94a-c09d-4a8d-a278-1109685aa98a.jpg)
![save_chem](https://user-images.githubusercontent.com/111716358/202489240-bd2e2070-6eab-43ba-b0e9-b4090185613c.jpg)

it can also visualize the result of well.

![Qa Va](https://user-images.githubusercontent.com/111716358/202493304-55de5209-766b-4b97-8e5c-09cf4273c736.jpg)
![Qg Vg](https://user-images.githubusercontent.com/111716358/202493341-b242c72f-b502-4d30-a01b-21040f36592c.jpg)
