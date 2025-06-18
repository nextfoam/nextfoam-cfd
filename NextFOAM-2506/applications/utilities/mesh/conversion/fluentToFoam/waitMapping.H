int waitMapping(uint16_t port)
{
    int clientSocket;
    struct timeval timeout;
    int rc;

    // creating socket

#if defined(WIN32)
    WSADATA wsadata;
    if (WSAStartup(MAKEWORD(1,1), &wsadata) == SOCKET_ERROR)
    {
        printf("Error creating socket.");
        return -1;
    }
#endif

    clientSocket = socket(AF_INET, SOCK_STREAM, 0);

    // set timeout for connection

    timeout.tv_sec = 1;
    timeout.tv_usec = 0;
    setsockopt
    (
        clientSocket, 
        SOL_SOCKET, 
        SO_SNDTIMEO, 
        (char *)&timeout, 
        sizeof(timeout)
    );

    // specifying address
    sockaddr_in serverAddress;
    serverAddress.sin_family = AF_INET;
    serverAddress.sin_port = htons(port);
    serverAddress.sin_addr.s_addr = inet_addr("127.0.0.1");

    // sending connection request

    //cout << "Connecting..." << endl;
    rc = connect
    (
        clientSocket, 
        (struct sockaddr*)&serverAddress,
        sizeof(serverAddress)
    );

    if (rc < 0)
        return -1;

    // sending data

    // cout << "Sending..." << endl;
    const char* message = "Write mapping";
    rc = send(clientSocket, message, strlen(message), 0);
    if (rc < 0)
        return -1;

    // Waiting data

    // cout << "Receiving..." << endl;
    char buffer[1024] = { 0 };
    recv(clientSocket, buffer, sizeof(buffer), 0);
    if (rc < 0)
        return -1;
    // cout << "Message from client: " << buffer << endl;

    // closing socket
#if defined(WIN32)
    closesocket(clientSocket);
    WSACleanup();
#else
    close(clientSocket);
#endif

    return 0;
}
